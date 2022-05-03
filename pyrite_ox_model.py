import boudreau_2010 as bd
from esbmtk import Source, Signal, Connect, DataField, ExternalData
import numpy as np

import matplotlib.pyplot as plt
import pandas as pd
import os

run_time = "200 kyr"
# run_time = "40 kyr"
time_step = "0.1 yr"
# rain_ratio = 0.238548387096774
# alpha = 0.495161290322581

rain_ratio = 0.3
alpha = 0.6

sum_dic: float = (2153 * 1.76e16) + (1952 * 2.85e16) + (2291 * 1.29e18)
total_vol: float = (1.76e16 + 2.85e16 + 1.29e18) * 1000  # m^3 -> L
dic: float = sum_dic / total_vol
sum_ta: float = (2345 * 1.76e16) + (2288 * 2.85e16) + (2399 * 1.29e18)  # mmol
ta: float = sum_ta / total_vol

M1 = bd.initialize_esbmtk_model(rain_ratio, alpha, run_time, time_step, ta, dic)


M1.Carbon_Pulse = Source(name="Carbon_Pulse", species=M1.CO2)

# #-------------------------------SCENARIO 1------------------------------------
# # 2CaCO3 + H2SO4 → 2Ca2+ + SO42- + 2HCO3-
# # Contributes +4 mol DIC, +4 mol TA
#
# M1.PyOx = Signal(name="PyOx",  # name of signal
#                species=M1.CO2,  # species
#                filename="koelling_model_200kyr.csv", # spreadsheet data is PYRITE OXIDATION flux
#                scale=4, # 1 mol FeS2 contributes 4 mol DIC
#                )
#
# Connect(
#     source=M1.Fw.DIC,
#     sink=M1.L_b.DIC,
#     rate="0 mol/yr",
#     signal=M1.PyOx,  # list of processes
#     id="PyOx_DIC"
# )
#
# Connect(
#     source=M1.Fw.TA,
#     sink=M1.L_b.TA,
#     ctype="scale_with_flux",
#     ref_reservoirs=M1.C_DIC_2_DIC_PyOx_DIC.DIC_2_DIC_PyOx_DIC_F, # this is the flux name
#     scale=1, # 1 mol FeS2 -> 2 mol CO2, 1 mol FeS2 -> 4 mol TA (note: Signal was already scaled by 2)
#     id="PyOx_TA",  # set a new id, e..g. PyOx_TA
# )

# #-------------------------------SCENARIO 2------------------------------------
# # H2SO4 + CaCO3 → SO42- + Ca2+ + CO2 + H2O
# # Contributes +2 mol DIC
# #
# # NOTE: Comment out the appropriate lines (do not need to plot TA) if using the Matplotlib plots!
# # See plot 'DIC/TA FLUX DATA'
#
# M1.PyOx = Signal(name="PyOx",  # name of signal
#                species=M1.CO2,  # species
#                filename="koelling_model_200kyr.csv", # spreadsheet data is PYRITE OXIDATION flux
#                scale=2, # 1 mol FeS2 contributes 2 mol DIC
#                )
#
# Connect(
#     source=M1.Fw.DIC,
#     sink=M1.L_b.DIC,
#     rate="0 mol/yr",
#     signal=M1.PyOx,  # list of processes
#     id="PyOx_DIC"
# )

#-------------------------------SCENARIO 3------------------------------------
# FeS2 + 14Fe3+ + 8H2O → 15Fe2+ + 2SO42- + 16H+
# Contributes -4 mol TA

# NOTE: Comment out the appropriate lines (do not need to plot DIC) if using the Matplotlib plots!
# See plot 'DIC/TA FLUX DATA'

M1.PyOx = Signal(name="PyOx",  # name of signal
               species=M1.CO2,  # species
               filename="koelling_model_200kyr.csv", # spreadsheet data is PYRITE OXIDATION flux
               scale=-4, # 1 mol FeS2 decreases TA by 4 mol
               )

Connect(
    source=M1.Fw.TA,
    sink=M1.L_b.TA,
    rate="0 mol/yr",
    signal=M1.PyOx,  # list of processes
    id="PyOx_TA"
)

M1.read_state()
M1.run(solver="numba")
#M1.save_state()
#M1.sub_sample_data()
M1.save_data()

#-------------------------------------------------------------------------------
fig, axs = plt.subplots(2, 3, figsize=(13,7))
time_x = M1.time[::-1] / 1000 # model run time for x axis

#-----------------------------DIC DATA------------------------------------------
axs[0,0].set_xlim(200, 0)
axs[0,0].set_ylim(1.8, 2.4)

axs[0,0].plot(time_x, M1.D_b.DIC.c * 1000, color="C0", label="Deep box")
axs[0,0].plot(time_x, M1.H_b.DIC.c * 1000, color="C1", label="High box")
axs[0,0].plot(time_x, M1.L_b.DIC.c  * 1000, color="C2", label="Low box")

# axs[0,0].set_title("DIC")
axs[0,0].set_xlabel("Time [ka]")
axs[0,0].set_ylabel("DIC [mmol/L]")
axs[0,0].legend(loc=5, prop={'size': 8})

# Marine isotope stages
axs[0,0].axvspan(0, 14, alpha=0.2, color='k')
axs[0,0].axvspan(29, 57, alpha=0.2, color='k')
axs[0,0].axvspan(71, 130, alpha=0.2, color='k')
axs[0,0].axvspan(191, 243, alpha=0.2, color='k')

y = axs[0,0].get_ylim()[1]
axs[0,0].annotate('MIS', (0, y), (8, -10), textcoords='offset points')
axs[0,0].annotate('1', (7.5, y), (0, -10), textcoords='offset points')
axs[0,0].annotate('2', (21.5, y), (0, -10), textcoords='offset points')
axs[0,0].annotate('3', (43, y), (0, -10), textcoords='offset points')
axs[0,0].annotate('4', (67, y), (0, -10), textcoords='offset points')
axs[0,0].annotate('5', (100.5, y), (0, -10), textcoords='offset points')
axs[0,0].annotate('6', (160.5, y), (0, -10), textcoords='offset points')
axs[0,0].annotate('7', (198, y), (0, -10), textcoords='offset points')

axs[0,0].annotate('(a)', (200, y), (-50, -10), textcoords='offset points', weight='bold')

#------------------------------TA DATA------------------------------------------
axs[0,1].set_xlim(200, 0)

axs[0,1].plot(time_x, M1.D_b.TA.c * 1000, color="C0", label="Deep box")
axs[0,1].plot(time_x, M1.H_b.TA.c * 1000, color="C1", label="High box")
axs[0,1].plot(time_x, M1.L_b.TA.c * 1000, color="C2", label="Low box")

# axs[0,1].set_title("TA")
axs[0,1].set_xlabel("Time [ka]")
axs[0,1].set_ylabel("Alkalinity [mmol/L]")
axs[0,1].legend(loc=5, prop={'size': 8})

y = axs[0,0].get_ylim()[1]
axs[0,1].axvspan(0, 14, alpha=0.2, color='k')
axs[0,1].axvspan(29, 57, alpha=0.2, color='k')
axs[0,1].axvspan(71, 130, alpha=0.2, color='k')
axs[0,1].axvspan(191, 243, alpha=0.2, color='k')

y = axs[0,1].get_ylim()[1]
axs[0,1].annotate('MIS', (0, y), (8, -10), textcoords='offset points')
axs[0,1].annotate('1', (7.5, y), (0, -10), textcoords='offset points')
axs[0,1].annotate('2', (21.5, y), (0, -10), textcoords='offset points')
axs[0,1].annotate('3', (43, y), (0, -10), textcoords='offset points')
axs[0,1].annotate('4', (67, y), (0, -10), textcoords='offset points')
axs[0,1].annotate('5', (100.5, y), (0, -10), textcoords='offset points')
axs[0,1].annotate('6', (160.5, y), (0, -10), textcoords='offset points')
axs[0,1].annotate('7', (198, y), (0, -10), textcoords='offset points')

axs[0,1].annotate('(b)', (200, y), (-50, -10), textcoords='offset points', weight='bold')

#-------------------------------pH DATA-----------------------------------------
axs[0,2].set_xlim(200, 0)

axs[0,2].plot(time_x, -np.log10(M1.D_b.cs.H), color="C0", label="Deep box")

axs[0,2].plot(time_x, -np.log10(M1.H_b.cs.H), color="C1", label="High box")

axs[0,2].plot(time_x, -np.log10(M1.L_b.cs.H), color="C2", label="Low box")

# axs[0,2].set_title("pH")
axs[0,2].set_xlabel("Time [ka]")
axs[0,2].set_ylabel("pH")
axs[0,2].legend(loc=5, prop={'size': 8})

axs[0,2].axvspan(0, 14, alpha=0.2, color='k')
axs[0,2].axvspan(29, 57, alpha=0.2, color='k')
axs[0,2].axvspan(71, 130, alpha=0.2, color='k')
axs[0,2].axvspan(191, 243, alpha=0.2, color='k')

y = axs[0,2].get_ylim()[1]
axs[0,2].annotate('MIS', (0, y), (8, -10), textcoords='offset points')
axs[0,2].annotate('1', (7.5, y), (0, -10), textcoords='offset points')
axs[0,2].annotate('2', (21.5, y), (0, -10), textcoords='offset points')
axs[0,2].annotate('3', (43, y), (0, -10), textcoords='offset points')
axs[0,2].annotate('4', (67, y), (0, -10), textcoords='offset points')
axs[0,2].annotate('5', (100.5, y), (0, -10), textcoords='offset points')
axs[0,2].annotate('6', (160.5, y), (0, -10), textcoords='offset points')
axs[0,2].annotate('7', (198, y), (0, -10), textcoords='offset points')

axs[0,2].annotate('(c)', (200, y), (-50, -10), textcoords='offset points', weight='bold')

#------------------------------DEPTHS DATA--------------------------------------
axs[1,0].set_xlim(200, 0)
axs[1,0].set_ylim(-4700, -3200)

axs[1,0].plot(time_x, -M1.D_b.cs.zsat, color="C0", label="zsat")

axs[1,0].plot(time_x, -M1.D_b.cs.zcc, color="C1", label="zcc")

axs[1,0].plot(time_x, -M1.D_b.cs.zsnow, color="C2", label="zsnow")

# axs[1,0].set_title("Depths")
axs[1,0].set_xlabel("Time [ka]")
axs[1,0].set_ylabel("Depth [m]")
axs[1,0].legend(loc=5, prop={'size': 8})

axs[1,0].axvspan(0, 14, alpha=0.2, color='k')
axs[1,0].axvspan(29, 57, alpha=0.2, color='k')
axs[1,0].axvspan(71, 130, alpha=0.2, color='k')
axs[1,0].axvspan(191, 243, alpha=0.2, color='k')

y = axs[1,0].get_ylim()[1]
axs[1,0].annotate('MIS', (0, y), (8, -10), textcoords='offset points')
axs[1,0].annotate('1', (7.5, y), (0, -10), textcoords='offset points')
axs[1,0].annotate('2', (21.5, y), (0, -10), textcoords='offset points')
axs[1,0].annotate('3', (43, y), (0, -10), textcoords='offset points')
axs[1,0].annotate('4', (67, y), (0, -10), textcoords='offset points')
axs[1,0].annotate('5', (100.5, y), (0, -10), textcoords='offset points')
axs[1,0].annotate('6', (160.5, y), (0, -10), textcoords='offset points')
axs[1,0].annotate('7', (198, y), (0, -10), textcoords='offset points')

axs[1,0].annotate('(d)', (200, y), (-50, -10), textcoords='offset points', weight='bold')

#-------------------------- pCO2 DATA-------------------------------------------
axs[1,1].set_xlim(200, 0)

axs[1,1].plot(time_x, M1.CO2_At.c * 1e6, label="pCO2") #mol/L to ppm

# axs[1,1].set_title("pCO2 in the Atmosphere")
axs[1,1].set_xlabel("Time [ka]")
axs[1,1].set_ylabel("pCO$_{2}$ [ppm]")

axs[1,1].axvspan(0, 14, alpha=0.2, color='k')
axs[1,1].axvspan(29, 57, alpha=0.2, color='k')
axs[1,1].axvspan(71, 130, alpha=0.2, color='k')
axs[1,1].axvspan(191, 243, alpha=0.2, color='k')

y = axs[1,1].get_ylim()[1]
axs[1,1].annotate('MIS', (0, y), (8, -10), textcoords='offset points')
axs[1,1].annotate('1', (7.5, y), (0, -10), textcoords='offset points')
axs[1,1].annotate('2', (21.5, y), (0, -10), textcoords='offset points')
axs[1,1].annotate('3', (43, y), (0, -10), textcoords='offset points')
axs[1,1].annotate('4', (67, y), (0, -10), textcoords='offset points')
axs[1,1].annotate('5', (100.5, y), (0, -10), textcoords='offset points')
axs[1,1].annotate('6', (160.5, y), (0, -10), textcoords='offset points')
axs[1,1].annotate('7', (198, y), (0, -10), textcoords='offset points')

axs[1,1].annotate('(e)', (200, y), (-50, -10), textcoords='offset points', weight='bold')

# # -------------------------- DIC/TA FLUX DATA----------------------------------
# axs[1,2].set_xlim(200, 0)
# # comment below out if running scenario 3
# pyox_dic = M1.C_DIC_2_DIC_PyOx_DIC.DIC_2_DIC_PyOx_DIC_F.m
# axs[1,2].plot(time_x, pyox_dic / 1E12, color="C0", label="DIC") # mol/yr -> Tmol/yr
#
# # comment below out if running scenario 2
# pyox_ta = M1.C_TA_2_TA_PyOx_TA.TA_2_TA_PyOx_TA_F.m
# axs[1,2].plot(time_x, pyox_ta / 1E12, color="C1", label="TA") # mol/yr -> Tmol/yr
#
# #axs[1,2].set_title("Flux")
# axs[1,2].set_xlabel("Time [ka]")
# axs[1,2].set_ylabel("Flux rate [Tmol/yr]")
# axs[1,2].legend(loc=3, prop={'size': 8})
#
# axs[1,2].axvspan(0, 14, alpha=0.2, color='k')
# axs[1,2].axvspan(29, 57, alpha=0.2, color='k')
# axs[1,2].axvspan(71, 130, alpha=0.2, color='k')
# axs[1,2].axvspan(191, 243, alpha=0.2, color='k')
#
# y = axs[1,2].get_ylim()[1]
# axs[1,2].annotate('MIS', (0, y), (8, -10), textcoords='offset points')
# axs[1,2].annotate('1', (7.5, y), (0, -10), textcoords='offset points')
# axs[1,2].annotate('2', (21.5, y), (0, -10), textcoords='offset points')
# axs[1,2].annotate('3', (43, y), (0, -10), textcoords='offset points')
# axs[1,2].annotate('4', (67, y), (0, -10), textcoords='offset points')
# axs[1,2].annotate('5', (100.5, y), (0, -10), textcoords='offset points')
# axs[1,2].annotate('6', (160.5, y), (0, -10), textcoords='offset points')
# axs[1,2].annotate('7', (198, y), (0, -10), textcoords='offset points')
#
# axs[1,2].annotate('(f)', (200, y), (-50, -10), textcoords='offset points', weight='bold')

#-----------------------------Adjusting Subplots--------------------------------
plt.rcParams["font.size"] = 12
fig.tight_layout()
plt.show()

# fig.savefig("pco2_plots/s1_4dic_4ta.png")
# fig.savefig("pco2_plots/s2_2dic_0ta.png")
# fig.savefig("pco2_plots/s3_0dic_4ta.png")


#-----------------PLOTTING USING DATAFIELDS------------------------------

# print(f"CO2 to H_b = {max(M1.C_H_b_2_CO2_At.CO2_F.m)/1E12:.2f} Tmol/yr")
# print(f"CO2 from L_b = {max(M1.C_L_b_2_CO2_At.CO2_F.m)/1E12:.2f} Tmol/yr")
# # print(f"Ex flux {max(M1.ASGE_DIC_2_CO2_At_EX.DIC_2_CO2_At_EX.m)/(1E12 * M1.dt):.2f}")
#
# DataField(
#     name="pH",
#     associated_with=M1.L_b.cs,
#     y1_data=[-np.log10(M1.L_b.cs.H), -np.log10(M1.H_b.cs.H), -np.log10(M1.D_b.cs.H)],
#     y1_legend=["Low latitude", "High latitude", "Deep box"],
#     y1_label="pH",
# )
#
# DataField(
#     name="DIC",
#     associated_with=M1.D_b.DIC,
#     y1_data=[M1.L_b.DIC.c * 1000, M1.H_b.DIC.c * 1000, M1.D_b.DIC.c * 1000],
#     y1_legend=["Low latitude", "High latitude", "Deep box"],
#     y1_label="DIC [mmol/L]",
# )
#
# DataField(
#     name="TA",
#     associated_with=M1.D_b.TA,
#     y1_data=[M1.L_b.TA.c * 1000, M1.H_b.TA.c * 1000, M1.D_b.TA.c * 1000],
#     y1_legend=["Low latitude", "High latitude", "Deep box"],
#     y1_label="TA [mmol/L]",
# )
#
# DataField(
#     name="CO3",
#     associated_with=M1.L_b.cs,
#     y1_data=[M1.L_b.cs.CO3 * 1000, M1.H_b.cs.CO3 * 1000, M1.D_b.cs.CO3 * 1000],
#     y1_legend=["Low latitude", "High latitude", "Deep box"],
#     y1_label="CO3 (mmol/L)",
# )
#
# DataField(
#     name="CO2aq",
#     associated_with=M1.L_b.cs,
#     y1_data=[M1.L_b.cs.CO2aq * 1000, M1.H_b.cs.CO2aq * 1000, M1.D_b.cs.CO2aq * 1000],
#     y1_legend=["Low latitude", "High latitude", "Deep box"],
#     y1_label="CO2 aq (mmol/L)",
# )
#
# DataField(
#     name="Depths",
#     associated_with=M1.D_b.cs,
#     y1_data=[-M1.D_b.cs.zsat, -M1.D_b.cs.zcc, -M1.D_b.cs.zsnow],
#     y1_legend=["zsat", "zcc", "zsnow"],
#     y1_label="Depth (m)",
# )
#
# dic_caw = M1.C_DIC_2_DIC_Ca_W.DIC_2_DIC_Ca_W_F.m
# dic_siw = M1.C_CO2_At_2_DIC_Si_W.CO2_At_2_DIC_Si_W_F.m
#
# ta_caw = M1.C_TA_2_TA_wca_ta.TA_2_TA_wca_ta_F.m
# ta_siw = M1.C_TA_2_TA_wsi_ta.TA_2_TA_wsi_ta_F.m
#
# dic_w = dic_caw + dic_siw
# ta_w = ta_caw + ta_siw
#
# # fw_dic = (M1.Caw.DIC_2_DIC_F.m + M1.Siw.CO2_At_2_DIC_F.m + M1.volcanic.DIC_2_CO2_At_F.m)/1e12
#
# DataField(
#     name="fluxes",
#     associated_with=M1.D_b.cs,
#     y1_data=[dic_w, ta_w],
#     y1_legend=["DIC_w", "TA_w"],
#     y1_label="Flux (mol/yr)",
# )
# # print(f"len of CaCO3 weathering DIC ={len(M1.Caw.DIC_2_DIC_F.m)}")
# # print(f"len of SiCO3 weathering DIC ={len(M1.Siw.CO2_At_2_DIC_F.m)}")
# # print(f"len of M1.D_b.cs.Fburial ={len(M1.D_b.cs.Fburial)}")
#
#
# # DataField(
# #     name="flux_D",
# #     associated_with=M1.D_b.cs,
# #     y1_data=((M1.D_b.cs.Fburial / 1e12) - (M1.CG_Fw2L_b.DIC.weathering_DIC_F.m / 1e12)),
# #     y1_legend="Fburial-Falk",
# #     y1_label="Flux (Tmol/yr)",
# # )
#
# DataField(
#     name="atm",
#     associated_with=M1.D_b.cs,
#     y1_data=M1.CO2_At.c * 1e6,
#     y1_legend="CO2",
#     y1_label="ppm",
# )
#
# # # Reading our digitized files
# # ExternalData(
# #     name="d_DIC",
# #     filename="D_DIC.csv",
# #     legend="D",
# #     # offset = "1838 yrs",
# #     reservoir=M1.D_b.DIC,
# # )
#
# # DataField(
# #     name="Carbon_pulse",
# #     associated_with=M1.D_b.cs,
# #     y1_data=M1.CP.data.m,
# #     y1_legend="CO2 input",
# #     y1_label="C mol/yr",
# # )
#
# # DataField(
# #     name="dd",
# #     associated_with=M1.D_b.DIC,
# #     y1_data=M1.d_DIC.y,
# #     x1_data=M1.d_DIC.x,
# #     y1_legend="CO2",
# #     y1_label="ppm",
# # )
#
# M1.plot([pH, DIC, TA, CO3, Depths, atm, fluxes])
# #M1.plot([Carbon_Pulse])
#
# diff = M1.C_H_b_2_CO2_At.CO2_F.m + M1.C_L_b_2_CO2_At.CO2_F.m
# print(f"CO2 uptake: {sum(diff[19900:20000]*M1.dt)*12/1E15:.2f}Pg")
#
# # CM.plot([dd])
# # CM.plot([pH, DIC, TA, CO3, Depths, fluxes, flux_D, atm])
#
#
# # print(f"co3 low lat box = {L_b.cs.CO3}")
# # print(f"co3 high lat box = {H_b.cs.CO3}")
# # print(f"co3 deep box = {D_b.cs.CO3}")
