#!/usr/bin/env python3
      
import os, sys
from netCDF4 import Dataset
from itertools import islice
# This script must be run in an environment with netCDF4 python libraries 
# active and the ability to build cesm cases.
# On Cheyenne, the best way to do this is from an interactive session.
# > qsub -X -I -l select=1:ncpus=36:mpiprocs=36 -l walltime=02:30:00 -q regular -A P93300642
# > conda activate my_env
# > ./build_ppe_CAM_cases.py

# Edit below to set up your cases

#cesmroot  = '/glade/work/addisus/sandboxes/cam6_4_127/'
cesmroot = '/glade/work/hannay/cesm_tags/cesm3_0_alpha07g'
basecasename = "wave1"
baseroot = os.path.join("/glade/work/qingyuany/simulations/", "cam7_correct_v", basecasename)
res = "ne30pg3_ne30pg3_mg17"
compset = "FHISTC_LTso" #"FHIST"

paramfile = "/glade/work/addisus/scam_scripts/PPE_250_ensemble/parameter38_100.nc" #"est_paras_scale20_original_scale.nc" #"parameter38_100.nc" 
ensemble_startval = "001" # The startval strings should be the same length, or else. 
basecase_startval = "000"
project = "P93300313"
wall_time = "12:00:00"
specify_user_num_sims = False # If True, use "user_num_sims" number of sims. This will
                             # use the first X values of each parameter, and will 
                             # probably crash if you specify more simulations than
                             # parameter values in your paramfile. 
                             # If False, automatically create the same number of cases
                             # as param values in the paramfile. This is the safest 
                             # option.
user_num_sims = 100

#command = "./xmlchange  -id CAM_CONFIG_OPTS -val ' -phys cam7'"
# Currently only number of years set (no months or days)
# Will be the same in every case

stop_yrs = 1
rest_yrs = 1

# REFCASE settings below are optional
# Will be the same in every case
# ref_case_get = False
# ref_case_name = "case.std"
# ref_case_path = "cesm2_init"

## CAM_CONFIG_OPTS
cam_conopts         = "-cosp"
#cam_phys_opt        = "-phys cam7"
#cam_phys_opt        = "-phys cam_dev"

## User_nl_cam text that will be cloned in every ensemble case
## Text in user_nl_cam will appear exactly as in this string specifier
# 'emulated', 'kk2000'
user_nl_string = """
cosp_passive=.true.

&camexp

    mfilt               =       0,       0,     20,      40,      12,       120,      1,   1
    nhtfrq              =       0,       0,    -24,      -3,       0,       -2,      0,  -8760
    ndens               =       2,       2,      2,       2,       2,       1,      2,   1
    interpolate_output  =  .true.,  .false., .true., .true., .false., .true.,  .true.
    interpolate_nlat    =     192,     192,    192,     192,     192,     192,   192
    interpolate_nlon    =     288,     288,    288,     288,     288,     288,   288 
    
    empty_htapes = .true.

    fincl1 = 'ACTNI', 'ACTNL', 'ACTREI', 'ACTREL', 'AODDUST', 'AODVIS', 'AODVISdn','BURDENBC',
    'BURDENDUST', 'BURDENPOM', 'BURDENSEASALT',
    'BURDENSO4', 'BURDENSOA', 'CAPE', 'CCN3', 'CDNUMC', 'CH4', 'CLDHGH', 'CLDICE', 'CLDLIQ', 'CLDLOW',
    'CLDMED', 'CLDTOT', 'CLOUD', 'CMFMC_DP',
    'CT_H2O', 'DCQ', 'DQCORE', 'DTCOND', 'DTCORE', 'DTV', 'EVAPPREC', 'EVAPSNOW', 'FCTI', 'FCTL', 'FICE', 'FLDS', 'FLNS', 'FLNSC', 'FLNT', 'FLNTC', 'FLUT',
    'FREQZM', 'FSDS', 'FSDSC', 'FSNS', 'FSNSC', 'FSNT', 'FSNTC', 'FSNTOA', 'ICEFRAC', 'LANDFRAC', 'LHFLX', 'LWCF', 'MPDICE', 'MPDLIQ', 'MPDQ', 'MPDT',
    'OCNFRAC', 'OMEGA', 'OMEGA500', 'PBLH', 'PHIS', 'PINT', 'PMID', 'PRECC', 'PRECL', 'PRECSC', 'PRECSL', 'PRECT', 'PS', 'PSL', 'PTEQ', 'PTTEND', 'Q',
    'QFLX', 'QRL', 'QRS', 'QTGW', 'RCMTEND_CLUBB', 'RELHUM', 'RVMTEND_CLUBB', 'SHFLX', 'SOLIN', 'SST',
    'STEND_CLUBB', 'SWCF', 'FLUTC',
    'T', 'TAUX', 'TAUY', 'TFIX', 'TGCLDIWP', 'TGCLDLWP', 'TMQ', 'TREFHT', 'TS', 'TTGW', 'U', 'U10',
    'UBOT', 'UTGWORO', 'UTGW_TOTAL',
    'V', 'VBOT', 'VTGWORO', 'VTGW_TOTAL', 'WPRTP_CLUBB', 'WPTHLP_CLUBB', 'Z3', 'ZMDQ', 'ZMDT', 'N2O',
    'CO2','CFC11','CFC12',
    'AODVISdn','AODDUSTdn','CCN3', 'CDNUMC', 'H2O', 'NUMICE', 'NUMLIQ','OMEGA500',
    'AQSO4_H2O2','AQSO4_O3', 'bc_a1', 'bc_a4', 'dst_a1', 'dst_a2', 'dst_a3', 'ncl_a1',
    'ncl_a1', 'ncl_a2', 'ncl_a3', 'pom_a1', 'pom_a4', 'so4_a1', 'so4_a2', 'so4_a3', 'soa_a2' ,
    'soa_a1', 'num_a1', 'num_a2', 'num_a3', 'num_a4',
    'bc_a1SFWET', 'bc_a4SFWET', 'dst_a1SFWET', 'dst_a2SFWET', 'dst_a3SFWET', 'ncl_a1SFWET',
    'ncl_a2SFWET', 'ncl_a3SFWET', 'pom_a1SFWET', 'pom_a4SFWET', 'so4_a1SFWET', 'so4_a2SFWET', 'so4_a3SFWET', 'soa_a1SFWET',
    'soa_a2SFWET', 'bc_c1SFWET', 'bc_c4SFWET', 'dst_c1SFWET', 'dst_c2SFWET', 'dst_c3SFWET', 'ncl_c1SFWET', 'ncl_c2SFWET',
    'ncl_c3SFWET', 'pom_c1SFWET', 'pom_c4SFWET', 'so4_c1SFWET', 'so4_c2SFWET', 'so4_c3SFWET', 'soa_c1SFWET', 'soa_c2SFWET',
    'bc_a1DDF', 'bc_a4DDF', 'dst_a1DDF', 'dst_a2DDF', 'dst_a3DDF', 'ncl_a1DDF', 'ncl_a2DDF', 'ncl_a3DDF',
    'pom_a1DDF', 'pom_a4DDF', 'so4_a1DDF', 'so4_a2DDF', 'so4_a3DDF', 'soa_a1DDF', 'soa_a2DDF',
    'so4_a1_CLXF', 'so4_a2_CLXF', 'SFbc_a4', 'SFpom_a4', 'SFso4_a1', 'SFso4_a2',
    'so4_a1_sfgaex1', 'so4_a2_sfgaex1', 'so4_a3_sfgaex1', 'soa_a1_sfgaex1', 'soa_a2_sfgaex1',
    'SFdst_a1','SFdst_a2', 'SFdst_a3', 'SFncl_a1', 'SFncl_a2', 'SFncl_a3',
    'num_a2_sfnnuc1', 'SFSO2', 'OCN_FLUX_DMS', 'SAD_SULFC', 'SAD_TROP', 'SAD_AERO', 'vIVT'

    fincl2 ='FISCCP1_COSP','CLDTOT_ISCCP','MEANCLDALB_ISCCP',
    'CLD_MISR','CLTMODIS','CLWMODIS','CLIMODIS','CLHMODIS','CLMMODIS','CLLMODIS',
    'LWPMODIS','IWPMODIS','CLMODIS'

    cosp_histfile_num=2

    fincl3 = 'PRECT', 'PRECC', 'FLUT', 'U850', 'U200', 'V850', 'V200', 'OMEGA500', 'TS', 'SST', 'PSL'
    fincl4 =  'PRECC','PRECL'
    fincl5 = 'Uzm','Vzm','Wzm','THzm', 'VTHzm','WTHzm','UVzm','UWzm'


    phys_grid_ctem_nfreq=-6
    phys_grid_ctem_zm_nbas=120
    phys_grid_ctem_za_nlat=90

    dust_emis_fact = 4.0
    seasalt_emis_scale = 1.0

    micro_mg_dcs = 500.D-6
    hetfrz_dust_scalfac = 0.20
    cldfrc_dp1 =  0.1
    clubb_C8 = 4.5

    se_sponge_del4_lev = 1
    se_sponge_del4_nu_div_fac = 1
    se_sponge_del4_nu_fac = 1
    se_hypervis_subcycle = 1
    se_nu_div = 1E15
    se_statefreq = 144

    ext_frc_specifier = 'bc_a4 -> /glade/campaign/acom/MUSICA/emis/cmip7/ne30/CEDS-CMIP-2025-04-18/bc_a4-em-AIR-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18_gn_175001-202312.nc',
    'num_a4 -> /glade/campaign/acom/MUSICA/emis/cmip7/ne30/CEDS-CMIP-2025-04-18/num_bc_a4-em-AIR-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18_gn_175001-202312.nc',
    'SO2 -> /glade/campaign/acom/MUSICA/emis/cmip7/ne30/CEDS-CMIP-2025-04-18/SO2-em-AIR-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18_gn_175001-202312.nc',
    'so4_a1 -> /glade/campaign/acom/MUSICA/emis/cmip7/ne30/CEDS-CMIP-2025-04-18/so4_a1_ene_vertical-em-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18_gn_175001-202312.nc',
    'num_a1 -> /glade/campaign/acom/MUSICA/emis/cmip7/ne30/CEDS-CMIP-2025-04-18/num_so4_a1_ene_vertical-em-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18_gn_175001-202312.nc',
         'num_a1 -> /glade/campaign/cesm/cesmdata/inputdata/atm/cam/chem/emis/historical_ne30pg3/emissions-cmip6_num_a1_so4_contvolcano_vertical_850-5000_ne30pg3_c20200125.nc',
         'num_a2 -> /glade/campaign/cesm/cesmdata/inputdata/atm/cam/chem/emis/historical_ne30pg3/emissions-cmip6_num_a2_so4_contvolcano_vertical_850-5000_ne30pg3_c20200125.nc',
         'so4_a1 -> /glade/campaign/cesm/cesmdata/inputdata/atm/cam/chem/emis/historical_ne30pg3/emissions-cmip6_so4_a1_contvolcano_vertical_850-5000_ne30pg3_c20200125.nc',
         'so4_a2 -> /glade/campaign/cesm/cesmdata/inputdata/atm/cam/chem/emis/historical_ne30pg3/emissions-cmip6_so4_a2_contvolcano_vertical_850-5000_ne30pg3_c20200125.nc',
         'SO2 -> /glade/campaign/cesm/cesmdata/inputdata/atm/cam/chem/emis/historical_ne30pg3/emissions-cmip6_SO2_contvolcano_vertical_850-5000_ne30pg3_c20200125.nc',

srf_emis_specifier = 'bc_a4 -> /glade/campaign/acom/MUSICA/emis/cmip7/ne30/CEDS-CMIP-2025-04-18/bc_a4-em-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18_gn_175001-202312.nc',
         'bc_a4 -> /glade/campaign/acom/MUSICA/emis/cmip7/ne30/DRES-CMIP-BB4CMIP7-2-0_smoothed_20251015/bc_a4_smoothed_input4MIPs_emissions_CMIP_DRES-CMIP-BB4CMIP7-2-0_gn_175001-202112.nc',
         'pom_a4 -> /glade/campaign/acom/MUSICA/emis/cmip7/ne30/CEDS-CMIP-2025-04-18/pom_a4-em-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18_gn_175001-202312.nc',
         'pom_a4 -> /glade/campaign/acom/MUSICA/emis/cmip7/ne30/DRES-CMIP-BB4CMIP7-2-0_smoothed_20251015/pom_a4_smoothed_input4MIPs_emissions_CMIP_DRES-CMIP-BB4CMIP7-2-0_gn_175001-202112.nc',
         'num_a4 -> /glade/campaign/acom/MUSICA/emis/cmip7/ne30/CEDS-CMIP-2025-04-18/num_bc_a4-em-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18_gn_175001-202312.nc',
         'num_a4 -> /glade/campaign/acom/MUSICA/emis/cmip7/ne30/DRES-CMIP-BB4CMIP7-2-0_smoothed_20251015/num_bc_a4_smoothed_input4MIPs_emissions_CMIP_DRES-CMIP-BB4CMIP7-2-0_gn_175001-202112.nc',
         'num_a4 -> /glade/campaign/acom/MUSICA/emis/cmip7/ne30/CEDS-CMIP-2025-04-18/num_pom_a4-em-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18_gn_175001-202312.nc',
         'num_a4 -> /glade/campaign/acom/MUSICA/emis/cmip7/ne30/DRES-CMIP-BB4CMIP7-2-0_smoothed_20251015/num_pom_a4_smoothed_input4MIPs_emissions_CMIP_DRES-CMIP-BB4CMIP7-2-0_gn_175001-202112.nc',
         'num_a1-> /glade/campaign/acom/MUSICA/emis/cmip7/ne30/CEDS-CMIP-2025-04-18/num_so4_a1_ag-em-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18_gn_175001-202312.nc',
         'num_a1-> /glade/campaign/acom/MUSICA/emis/cmip7/ne30/CEDS-CMIP-2025-04-18/num_so4_a1_ship_slv-em-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18_gn_175001-202312.nc',
         'num_a2-> /glade/campaign/acom/MUSICA/emis/cmip7/ne30/CEDS-CMIP-2025-04-18/num_so4_a2_res_trs-em-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18_gn_175001-202312.nc',
         'so4_a1-> /glade/campaign/acom/MUSICA/emis/cmip7/ne30/CEDS-CMIP-2025-04-18/so4_a1_ag_ship_slv-em-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18_gn_175001-202312.nc',
         'so4_a2-> /glade/campaign/acom/MUSICA/emis/cmip7/ne30/CEDS-CMIP-2025-04-18/so4_a2_res_trs-em-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18_gn_175001-202312.nc',
         'so4_a1-> /glade/campaign/acom/MUSICA/emis/cmip7/ne30/DRES-CMIP-BB4CMIP7-2-0_smoothed_20251015/SO4_smoothed_input4MIPs_emissions_CMIP_DRES-CMIP-BB4CMIP7-2-0_gn_175001-202112.nc',
         'num_a1-> /glade/campaign/acom/MUSICA/emis/cmip7/ne30/DRES-CMIP-BB4CMIP7-2-0_smoothed_20251015/num_SO4_smoothed_input4MIPs_emissions_CMIP_DRES-CMIP-BB4CMIP7-2-0_gn_175001-202112.nc',
         'SO2 -> /glade/campaign/acom/MUSICA/emis/cmip7/ne30/CEDS-CMIP-2025-04-18/SO2-em-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18_gn_175001-202312.nc',
         'SO2 -> /glade/campaign/acom/MUSICA/emis/cmip7/ne30/DRES-CMIP-BB4CMIP7-2-0_smoothed_20251015/SO2_smoothed_input4MIPs_emissions_CMIP_DRES-CMIP-BB4CMIP7-2-0_gn_175001-202112.nc',
         'SOAE -> 0.5954D0*/glade/campaign/acom/MUSICA/emis/cmip7/ne30/DRES-CMIP-BB4CMIP7-2-0_smoothed_20251015/ISOP_smoothed_input4MIPs_emissions_CMIP_DRES-CMIP-BB4CMIP7-2-0_gn_175001-202112.nc',
         'SOAE -> 2.5592D0*/glade/campaign/acom/MUSICA/emis/cmip7/ne30/DRES-CMIP-BB4CMIP7-2-0_smoothed_20251015/BENZENE_smoothed_input4MIPs_emissions_CMIP_DRES-CMIP-BB4CMIP7-2-0_gn_175001-202112.nc',
         'SOAE -> 2.5592D0*/glade/campaign/acom/MUSICA/emis/cmip7/ne30/CEDS-CMIP-2025-04-18/BENZENE-VOC-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18-supplemental_gn_175001-202312.nc',
         'SOAE -> 8.5371D0*/glade/campaign/acom/MUSICA/emis/cmip7/ne30/DRES-CMIP-BB4CMIP7-2-0_smoothed_20251015/IVOC_smoothed_input4MIPs_emissions_CMIP_DRES-CMIP-BB4CMIP7-2-0_gn_175001-202112.nc',
         'SOAE -> 8.5371D0*/glade/campaign/acom/MUSICA/emis/cmip7/ne30/CEDS-CMIP-2025-04-18/IVOC-em-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18_gn_175001-202312.nc',
         'SOAE -> 16.650D0*/glade/campaign/acom/MUSICA/emis/cmip7/ne30/DRES-CMIP-BB4CMIP7-2-0_smoothed_20251015/SVOC_smoothed_input4MIPs_emissions_CMIP_DRES-CMIP-BB4CMIP7-2-0_gn_175001-202112.nc',
         'SOAE -> 16.650D0*/glade/campaign/acom/MUSICA/emis/cmip7/ne30/CEDS-CMIP-2025-04-18/SVOC-em-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18_gn_175001-202312.nc',
         'SOAE -> 8.2367D0*/glade/campaign/acom/MUSICA/emis/cmip7/ne30/DRES-CMIP-BB4CMIP7-2-0_smoothed_20251015/TOLUENE_smoothed_input4MIPs_emissions_CMIP_DRES-CMIP-BB4CMIP7-2-0_gn_175001-202112.nc',
         'SOAE -> 8.2367D0*/glade/campaign/acom/MUSICA/emis/cmip7/ne30/CEDS-CMIP-2025-04-18/TOLUENE-VOC-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18-supplemental_gn_175001-202312.nc',
         'SOAE -> 6.5013D0*/glade/campaign/acom/MUSICA/emis/cmip7/ne30/DRES-CMIP-BB4CMIP7-2-0_smoothed_20251015/XYLENES_smoothed_input4MIPs_emissions_CMIP_DRES-CMIP-BB4CMIP7-2-0_gn_175001-202112.nc',
         'SOAE -> 6.5013D0*/glade/campaign/acom/MUSICA/emis/cmip7/ne30/CEDS-CMIP-2025-04-18/XYLENES-VOC-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18-supplemental_gn_175001-202312.nc',
         'SOAE -> 5.1004D0*/glade/campaign/acom/MUSICA/emis/cmip7/ne30/DRES-CMIP-BB4CMIP7-2-0_smoothed_20251015/MTERP_smoothed_input4MIPs_emissions_CMIP_DRES-CMIP-BB4CMIP7-2-0_gn_175001-202112.nc'

/


"""
## No /EOF needed for python

## Should not be a need to edit below this line ##
'''
if cesmroot is None:
    raise SystemExit("ERROR: CESM_ROOT must be defined in environment")
_LIBDIR = os.path.join(cesmroot,"cime","scripts","Tools")
sys.path.append(_LIBDIR)
_LIBDIR = os.path.join(cesmroot,"cime","scripts","lib")
sys.path.append(_LIBDIR)

'''
# For CAM_dev_pumas, the location for Tools is in $cesmroot/cime/CIME
if cesmroot is None:
    raise SystemExit("ERROR: CESM_ROOT must be defined in environment")
_LIBDIR = os.path.join(cesmroot,"cime", "CIME","Tools")
sys.path.append(_LIBDIR)
_LIBDIR = os.path.join(cesmroot,"cime")
sys.path.append(_LIBDIR)


print(_LIBDIR)

import datetime, glob, shutil
import CIME.build as build
import numpy as np
from standard_script_setup import *
from CIME.case             import Case
from CIME.utils            import safe_copy
from argparse              import RawTextHelpFormatter
from CIME.locked_files          import lock_file, unlock_file

def per_run_case_updates(case, ensemble_str, nint, paramdict, ensemble_idx):
    print(">>>>> BUILDING CLONE CASE...")
    caseroot = case.get_value("CASEROOT")
    basecasename = os.path.basename(caseroot)[:-nint]

    unlock_file("env_case.xml",caseroot=caseroot)
    casename = basecasename+ensemble_str
    case.set_value("CASE",casename)
    rundir = case.get_value("RUNDIR")
    rundir = rundir[:-nint]+ensemble_str
    case.set_value("RUNDIR",rundir)
    case.flush()
    lock_file("env_case.xml",caseroot=caseroot)
    print("...Casename is {}".format(casename))
    print("...Caseroot is {}".format(caseroot))
    print("...Rundir is {}".format(rundir))

    # Add user_nl updates for each run                                                        

    paramLines = []
#    ens_idx = int(ensemble_str)-int(ensemble_startval)
    ens_idx = int(ensemble_idx)-int(ensemble_startval)
    print('ensemble_index')
    print(ens_idx)
    for var in paramdict.keys():
        paramLines.append("{} = {}\n".format(var,paramdict[var][ens_idx]))

    usernlfile = os.path.join(caseroot,"user_nl_cam")
    print("...Writing to user_nl file: "+usernlfile)
    file1 = open(usernlfile, "a")
    file1.writelines(paramLines)
    file1.close()

    print(">> Clone {} case_setup".format(ensemble_str))
    case.case_setup()
    print(">> Clone {} create_namelists".format(ensemble_str))
    case.create_namelists()
    #Addisu comments the following two lines
    print(">> Clone {} submit".format(ensemble_str))
    case.submit()

#nuopc #moab #mct
#'cheyenne'
def build_base_case(baseroot, basecasename, res, compset, overwrite):
    print(">>>>> BUILDING BASE CASE...")
    caseroot = os.path.join(baseroot,basecasename+'.'+basecase_startval)
    if overwrite and os.path.isdir(caseroot):
        shutil.rmtree(caseroot)
    with Case(caseroot, read_only=False) as case:
        if not os.path.isdir(caseroot):
            case.create(os.path.basename(caseroot), cesmroot, compset, res,
                        machine_name="derecho", driver="nuopc", 
                        run_unsupported=True, answer="r",walltime=wall_time, 
                        project=project)
            # make sure that changing the casename will not affect these variables                                           
            case.set_value("EXEROOT",case.get_value("EXEROOT", resolved=True))
            case.set_value("RUNDIR",case.get_value("RUNDIR",resolved=True)+".00")

            case.set_value("RUN_TYPE","startup")
            #case.set_value("GET_REFCASE",ref_case_get)
            #case.set_value("RUN_REFCASE",ref_case_name)
            #case.set_value("RUN_REFDIR",ref_case_path)
            case.set_value("STOP_OPTION","nyears")
            case.set_value("STOP_N",stop_yrs)
            case.set_value("REST_OPTION","nyears")
            case.set_value("REST_N",rest_yrs)
            case.set_value("RUN_STARTDATE", "2004-01-01")
            #case.set_value( "CAM_CONFIG_OPTS",
            #               case.get_value("CAM_CONFIG_OPTS",resolved=True) + ' ' + cam_phys_opt)
            case.set_value("CAM_CONFIG_OPTS", 
                           case.get_value("CAM_CONFIG_OPTS",resolved=True)+' '+cam_conopts)


        rundir = case.get_value("RUNDIR")
        caseroot = case.get_value("CASEROOT")
        
        print(">> base case_setup...")
        case.case_setup()
        
        print(">> base case write user_nl_cam...")
        usernlfile = os.path.join(caseroot,"user_nl_cam")
        funl = open(usernlfile, "a")
        funl.write(user_nl_string)
        funl.close()
        
        print(">> base case_build...")
        build.case_build(caseroot, case=case)

        return caseroot

def clone_base_case(caseroot, ensemble, overwrite, paramdict, ensemble_num):
    print(">>>>> CLONING BASE CASE...")
    print(ensemble)
    startval = ensemble_startval
    nint = len(ensemble_startval)
    cloneroot = caseroot
    for i in range(int(startval), int(startval)+ensemble):
        ensemble_idx = '{{0:0{0:d}d}}'.format(nint).format(i)
        member_string = ensemble_num[i-1]
        print(member_string)
        member_string = str(member_string)
        #member_string = np.array(member_string, dtype=str)
        print("member_string="+member_string)
        if ensemble > 1:
            caseroot = caseroot[:-nint] + member_string
        if overwrite and os.path.isdir(caseroot):
            shutil.rmtree(caseroot)
        if not os.path.isdir(caseroot):
            with Case(cloneroot, read_only=False) as clone:
                clone.create_clone(caseroot, keepexe=True)
        with Case(caseroot, read_only=False) as case:
            per_run_case_updates(case, member_string, nint, paramdict, ensemble_idx)


def take(n, iterable):
    "Return first n items of the iterable as a list"
    return list(islice(iterable, n))


def _main_func(description):

    print ("Starting SCAM PPE case creation, building, and submission script")
    print ("Base case name is {}".format(basecasename))
    print ("Parameter file is "+paramfile)

    overwrite = True

    # read in NetCDF parameter file
    inptrs = Dataset(paramfile,'r')
    print ("Variables in paramfile:")
    print (inptrs.variables.keys())
    print ("Dimensions in paramfile:")
    print (inptrs.dimensions.keys())
    num_sims = inptrs.dimensions['nmb_sim'].size
    num_vars = len(inptrs.variables.keys())-1
    ensemble_num = inptrs['Sample_nmb']
    print('ensemble_num')
    print(ensemble_num)
    #print(ensemble_num[0])
    if specify_user_num_sims:
        num_sims = user_num_sims

    print ("Number of sims = {}".format(num_sims))
    print ("Number of params = {}".format(num_vars))


    # Save a pointer to the netcdf variables
    paramdict = inptrs.variables
    del paramdict['Sample_nmb']

    # Create and build the base case that all PPE cases are cloned from
    caseroot = build_base_case(baseroot, basecasename, res,
                            compset, overwrite)

    # Pass in a dictionary with all of the parameters and their values
    # for each PPE simulation 
    # This code clones the base case, using the same build as the base case,
    # Adds the namelist parameters from the paramfile to the user_nl_cam 
    # file of each new case, does a case.set up, builds namelists, and 
    # submits the runs.
    # Addisu comment this
    clone_base_case(caseroot, num_sims, overwrite, paramdict, ensemble_num)

    inptrs.close()

if __name__ == "__main__":
    _main_func(__doc__)

