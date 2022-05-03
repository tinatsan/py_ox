# Import this code into a model setup script. See e.g., readsignal2.py
def initialize_esbmtk_model(rain_ratio, alpha, run_time, time_step, ta=0, dic=0):
    from esbmtk import (
        Model,
        Q_,
        AirSeaExchange,
        GasReservoir,
        DataField,
        create_bulk_connections,
        create_reservoirs,
        build_ct_dict,
        add_carbonate_system_1,
        add_carbonate_system_2,
        ExternalData,
        Connect,
    )

    from esbmtk import phc
    import numpy as np
    import logging

    # print(f"rain = {rain_ratio}, alpha = {alpha}")

    M1 = Model(
        name="M1",  # model name
        stop=run_time,  # end time of model
        timestep=time_step,  # time step
        element=[
            "Carbon",
            "Boron",
            "Hydrogen",
            "Phosphor",
            "Oxygen",
        ],
        mass_unit="mol",
        volume_unit="liter",
        concentration_unit="mol/kg",
        isotopes=False,
        # step_limit=5e5,
    )

    # -------------------- Set up box parameters -------------------------
    a_H = 0.15  # Area percentage High box
    a_L = 0.85  # Area percentage Low box
    a_D = 1.00  # Area percentage Deep ocean

    """ boxes are defined by area and depth interval here we use an ordered
    dictionary to define the box geometries. The next column is temperature
    in deg C, followed by pressure in bar
    the geometry is [upper depth datum, lower depth datum, area percentage]
    """
    bn: dict = {  # name: [[geometry], T, P]
        # High-Lat Box
        "H_b": {"g": [-200, -550, a_H], "T": 2, "P": 17.6, "S": 35},  # box height 350 m
        # Low-Lat Box
        "L_b": {"g": [-200, -300, a_L], "T": 21.5, "P": 5, "S": 35},  # box height 100 m
        # Deep Box
        "D_b": {
            "g": [-130, -6000, a_D],
            "T": 2,
            "P": 240,
            "S": 35,
        },  # box height 3852 m
        "Fw": {"ty": "Source", "sp": [M1.DIC, M1.TA]},
        # "Bw": {"ty": "Sink", "sp": [M1.DIC, M1.TA]},
    }

    # Calculating initial DIC and TA:
    Hv = 1.76e16
    Lv = 2.85e16
    Dv = 1.29e18

    total_vol: float = Hv + Lv + Dv

    if dic == 0:
        dic = 2.98e18 / (total_vol * 1)

    if ta == 0:
        ta = 3.12e18 / (total_vol * 1)

    print(f"Initial Model TA = {ta}, DIC = {dic} ta/dic = {ta/dic}")
    """ Set up initial conditions. Here we use shortcut and use the
    same conditions in each box. You can also specify the full list,
    see the output lic. This dict additionally specifies whether a
    species will require isotope calculations. If you need box
    specific initial conditions use the output of
    build_concentration_dicts as starting point, but it is usually
    faster to just run the model to steady state and use this as a
    starting condition """

    ic: dict = {
        # species: concentration, Isotopes, delta value
        M1.DIC: [Q_(f"{dic} mmol/kg"), False, 0], # mmol/kg
        M1.TA: [Q_(f"{ta} mmol/kg"), False, 0], # mmol/kg
    }
    # initialize reservoirs
    create_reservoirs(bn, ic, M1)

    # Ocean mixing  Conveyor belt
    thc = Q_("25 Sv")

    # Specify upwelling connections
    thx = {
        "L_b_to_H_b@thc": thc * M1.L_b.swc.density/1000,
        # Downwelling
        "H_b_to_D_b@thc": thc * M1.H_b.swc.density/1000,
        # Upwelling
        "D_b_to_L_b@thc": thc * M1.D_b.swc.density/1000,
    }
    sl: list = list(ic.keys())  # get species list
    d1 = {"ty": "scale_with_concentration", "sp": sl}
    ct = build_ct_dict(thx, d1)

    create_bulk_connections(ct, M1, mt="1:1")

    # mixing fluxes
    mix = {
        # surface/intermediate water mixing
        "H_b_to_D_b@mix_down": Q_("30 Sv") * M1.H_b.swc.density/1000,
        "D_b_to_H_b@mix_up": Q_("30 Sv") * M1.D_b.swc.density/1000,
    }
    ct = build_ct_dict(mix, d1)
    create_bulk_connections(ct, M1, mt="1:1")

    """ Organic matter flux P - 200 Tm/yr POM = Particulate Organic
    matter. Unlike mix, this is really just a name to keep things
    organized Since this model uses a fixed rate we can declare this
    flux with the rate keyword. POM only affects DIC, so we write """

    OM_export = Q_("200 Tmol/yr")
    CaCO3_export = OM_export * rain_ratio

    # Fluxes going into deep box
    ct = {  # organic matter
        "L_b_to_D_b@POM": {
            "sp": M1.DIC,
            "ty": "Regular",
            "ra": OM_export,
        },  # DIC from CaCO3
        "L_b_to_D_b@PIC_DIC": {
            "sp": M1.DIC,
            "ty": "Regular",
            "ra": CaCO3_export,
            "bp": "sink",
        },
        # # TA from CaCO3
        "L_b_to_D_b@PIC_TA": {
            "sp": M1.TA,
            "ty": "Regular",
            "ra": CaCO3_export * 2,
            "bp": "sink",
        },
    }
    create_bulk_connections(ct, M1, mt="1:1")
    M1.ct = ct

    # ----- Carbonate System- virtual reservoir definitions -------

    surface_boxes: list = [M1.L_b, M1.H_b]
    deep_boxes: list = [M1.D_b]
    print(M1.flux_summary())
    ef = M1.flux_summary(filter_by="PIC_DIC", return_list=True)
    print("ef is")
    print(ef)

    add_carbonate_system_1(surface_boxes)

    add_carbonate_system_2(
        rgs=deep_boxes,  # list of reservoir groups
        carbonate_export_fluxes=ef[1],  # list of fluxes
        zsat_min=-200,  # zsat_max
        z0=-200,
        alpha=alpha,
    )

    # -------------------- Atmosphere -------------------------
    GasReservoir(
        name="CO2_At",
        species=M1.CO2,
        reservoir_mass="1.833E20 mol",
        species_ppm="270 ppm",
    )

    AirSeaExchange(
        gas_reservoir=M1.CO2_At,  # Reservoir
        liquid_reservoir=M1.H_b.DIC,  # ReservoirGroup
        species=M1.CO2,
        ref_species=M1.H_b.cs.CO2aq,
        solubility=M1.H_b.swc.SA_co2 / 1,  # float
        area=M1.H_b.area,
        piston_velocity="4.8 m/d",
        water_vapor_pressure=M1.H_b.swc.p_H2O,
        id="gex_hb",
    )

    AirSeaExchange(
        gas_reservoir=M1.CO2_At,  # Reservoir
        liquid_reservoir=M1.L_b.DIC,  # ReservoirGroup
        species=M1.CO2,
        ref_species=M1.L_b.cs.CO2aq,
        solubility=M1.L_b.swc.SA_co2 / 1,  # float
        area=M1.L_b.area,
        piston_velocity="4.8 m/d",
        water_vapor_pressure=M1.L_b.swc.p_H2O,
        id="gex_lb",
    )

    # -------------------- Sources and Signals -------------------------

    # volcanic flux
    Connect(
        source=M1.Fw.DIC,
        sink=M1.CO2_At,  # source of flux
        cype="Regular",
        name="volcanic",
        rate="5 / 17 * 1e13 mol/yr",
        id="volcanic_flux",
    )

    # CaCO3 weathering DIC Land to Ocean DIC
    Connect(
        source=M1.Fw.DIC,  # source of flux
        sink=M1.L_b.DIC,
        reservoir_ref=M1.CO2_At,
        ctype="weathering",
        id="Ca_W",
        scale=1,
        ex=0.4,
        pco2_0="280 ppm",
        rate="12 / 17 * 1e13 mol/yr",  # ~70% of total flux
    )

    # CaSiO3 weathering Atmosphere to Ocean DIC
    Connect(
        source=M1.CO2_At,
        sink=M1.L_b.DIC,  # source of flux
        reservoir_ref=M1.CO2_At,
        ctype="weathering",
        id="Si_W",
        scale=1,
        ex=0.2,
        pco2_0="280 ppm",
        rate="5 / 17 * 1e13 mol/yr",  # ~30% of total flux
    )

    # CaCO3 weathering Land to Ocean ALK
    wf = M1.flux_summary(filter_by="Ca_W", return_list=True)

    print(wf)
    print(M1.flux_summary())

    Connect(
        source=M1.Fw.TA,
        sink=M1.L_b.TA,
        ctype="scale_with_flux",
        ref_flux=wf[1],
        id="wca_ta",
        scale=2,
    )

    # CaSiO3 weathering Land to Ocean ALK
    wf = M1.flux_summary(filter_by="Si_W", return_list=True)

    Connect(
        source=M1.Fw.TA,
        sink=M1.L_b.TA,
        ctype="scale_with_flux",
        ref_flux=wf[1],
        id="wsi_ta",
        scale=2,
    )

    return M1


# this goes at the end of the file
if __name__ == "__main__":

    from esbmtk import DataField
    import sys

    # run_time = "200 kyr"
    run_time = "1000 yr"
    time_step = "0.1 yr"
    rain_ratio = 0.3
    alpha = 0.6
    # ta = 2.2241935483871
    # dic = 2.24193548387097
    sum_dic: float = (2153 * 1.76e16) + (1952 * 2.85e16) + (2291 * 1.29e18)
    total_vol: float = (1.76e16 + 2.85e16 + 1.29e18) * 1000  # m^3 -> L
    dic: float = sum_dic / total_vol
    sum_ta: float = (2345 * 1.76e16) + (2288 * 2.85e16) + (2399 * 1.29e18)  # mmol
    ta: float = sum_ta / total_vol

    M1 = initialize_esbmtk_model(rain_ratio, alpha, run_time, time_step, ta, dic)
    dic = M1.D_b.DIC.m[-2] + M1.H_b.DIC.m[-2] + M1.L_b.DIC.m[-2]
    ta = M1.D_b.TA.m[-2] + M1.H_b.TA.m[-2] + M1.L_b.TA.m[-2]
    print(f"ta/dic = {ta/dic}, ta = {ta:.2e}, dic = {dic:.2e}")

    M1.run(solver="numba")
    # M1.run()
    # M1.save_state()

    dic_caw = M1.C_DIC_2_DIC_Ca_W.DIC_2_DIC_Ca_W_F.fa[0]
    dic_siw = M1.C_CO2_At_2_DIC_Si_W.CO2_At_2_DIC_Si_W_F.fa[0]

    ta_caw = M1.C_TA_2_TA_wca_ta.TA_2_TA_wca_ta_F.fa[0]
    ta_siw = M1.C_TA_2_TA_wsi_ta.TA_2_TA_wsi_ta_F.fa[0]
    # ra_siw = M1.L_b.TA
    print(f"DIC from Carbonate Weathering = {dic_caw:.2e}")
    print(f"DIC from Silicate Weathering = {dic_siw:.2e}")
    print(f"Total DIC flux : {(dic_caw+dic_siw):.2e}\n")

    print(f"TA from Carbonate Weathering = {ta_caw:.2e}")
    print(f"TA from Silicate Weathering = {ta_siw:.2e}")
    print(f"Total TA flux : {(ta_caw+ta_siw):.2e}\n")

    # dic_w = (M1.Caw.DIC_2_DIC_F.m[3] + M1.Siw.CO2_At_2_DIC_F.m[3])/1e12
    # ta_w =  (M1.TA.TA_2_TA_wca_ta_F.m[3] + M1.wsi_ta.TA_2_TA_F.m[3])/1e12

    # print(f"Carbonate DIC: {M1.Caw.DIC_2_DIC_F.m[3]:.2e}, Carbonate TA: {M1.TA.TA_2_TA_wca_ta_F.m[3]:.2e}")
    # print(f"Silicate DIC: {M1.Siw.CO2_At_2_DIC_F.m[3]:.2e}, Silicate TA: {M1.wsi_ta.TA_2_TA_F.m[3]:.2e}")

    # print(f"DIC Weathering = {dic_w:.2f}, TA Weathering {ta_w:.2f}")

    # dic = M1.D_b.DIC.m + M1.H_b.DIC.m + M1.L_b.DIC.m
    # ta = M1.D_b.TA.m + M1.H_b.TA.m + M1.L_b.TA.m

    # print(f"CO2 to H_b = {M1.C_H_b_2_CO2_At.CO2_F.m[-2]/1E12:.2f}")
    # print(f"CO2 from L_b = {M1.C_L_b_2_CO2_At.CO2_F.m[-2]/1E12:.2f}")
    # # print(f"Ex flux {max(M1.ASGE_DIC_2_CO2_At_EX.DIC_2_CO2_At_EX.m)/(1E12 * M1.dt):.2f}")
    # print(f"ta/dic = {ta[-2]/dic[-2]}, ta = {ta[-2]:.2e}, dic = {dic[-2]:.2e}")M

    print(f"pCO2 = {M1.CO2_At.c[-2]*1E6:.0f} ppm")
    # # print(f"Final TA = {ta}, DIC = {dic} ta/dic = {ta/dic}")
    # # run_esbmtk_model(rain_ratio, alpha)
    DataField(
        name="Depths",
        associated_with=M1.D_b.cs,
        y1_data=[-M1.D_b.cs.zsat[3:-2], -M1.D_b.cs.zcc[3:-2], -M1.D_b.cs.zsnow[3:-2]],
        x1_data=[M1.time[3:-2], M1.time[3:-2], M1.time[3:-2]],
        y1_legend=["zsat", "zcc", "zsnow"],
        y1_label="Depth (m)",
    )
    # DataField(
    #     name="ratios",
    #     associated_with=M1.D_b.cs,
    #     y1_data=[ta / dic],
    #     y1_legend=["ta/dic"],
    #     y1_label="ta/dic",
    # )

    M1.plot([Depths])
