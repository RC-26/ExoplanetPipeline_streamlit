import numpy as np ; import math
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import sys ; import os ; import glob
import timeit
from PIL import Image
from ics import Calendar, Event ; import datetime
from timezonefinder import TimezoneFinder ; import pytz
import streamlit as st

import astropy
import astropy.units as u
import astropy.time as at
from astropy.io import fits
from astropy.wcs import WCS
from astropy.utils.data import get_pkg_data_filename
import astropy.utils.data as aud
from astropy.coordinates import SkyCoord, AltAz, EarthLocation, get_body
import astropy.coordinates as ac
import astropy.constants as constants

import astroquery.jplhorizons as jpl
from astroquery.simbad import Simbad
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive as NEA
from astroquery.vizier import Vizier

from astroplan import FixedTarget, Observer
import astroplan as ap
import astroplan.plots as app
from astroplan.plots import plot_sky

st.title('Observation Scheduling Tool')
st.markdown ('This pipeline is for helping users plan their observations by extracting data from the NASA Exoplanet Archive to predict the next transit events. As of 2025 December 22, this pipeline mainly focuses on the ground-based telescopes that are part of the SkyNet Robotic Telescope Network.')

# targets = 'TRAPPIST-1'

targets = st.text_input (label = 'Target planet: ', placeholder = 'Example: TRAPPIST-1 c')
t_cb    = st.checkbox ('Submit target planet')
if t_cb:
    start_date = st.datetime_input (label = 'Start date: ')
    sd_cb      = st.checkbox ('Submit start date')
    if sd_cb:
        end_date  = st.datetime_input (label = 'End date: ', min_value = start_date)
        ed_submit = st.checkbox ('Submit end date')
        if ed_submit:
            ####################################################################################################
            # SkyNet Observatories
            SN_obs = ['Cerro Tololo Inter-American Obs', 'Meckering Obs', 'Perth Obs', 'American Public University System Obs',
                      'Astronomical Obs of the Jagiellonian University', 'Fan Mountain Obs', 'Hampden-Sydney College Obs',
                      'Montana Learning Center', 'Morehead', 'Northern Skies Obs', 'Sleaford Obs']
            SN_long       = [-70.805, 116.989, 116.136, -77.863, 19.828, -78.694, -78.471, -105.53, -79.05, -72.166, -105.921]
            SN_lat        = [-30.168, -31.638, -32.007,  39.293, 50.054,  37.879,  37.238,  32.902, 35.914,  44.325,   52.085]
            SN_elev       = [   2286,     197,     386,     170,    318,     546,     164,    2225,    145,     384,      580]
            SN_aperture   = [     40,      16,      16,      24,     20,      24,      16,      16,     24,      16,       16]
            SN_telescopes = [
                ['CTIO-1.0m', 'Prompt1', 'Prompt2', 'Prompt3', 'Prompt5', 'Prompt6', 'Prompt7', 'Prompt8'],
                ['PROMPT-MO-1'],
                ['R-COP'],
                ['APUS-CDK24'],
                ['OAUJ-CDK500'],
                ['RRRT'],
                ['HSC'],
                ['MLC-RCOS16'],
                ['Morehead'],
                ['NSO-17-CDK'],
                ['PROMPT-USASK', 'PROMPT-USASK-2', 'USASK-14']
            ]
            
            def Timezone_Finder (lon, lat, elev):
                location = EarthLocation(lon =lon*u.deg, lat = lat*u.deg, height = elev*u.m)
                tf = TimezoneFinder()
                timezone_str = tf.certain_timezone_at(lng=location.lon.deg, lat=location.lat.deg)
                if timezone_str is None:
                    print("Could not determine the time zone for the given location.")
                else:
                    t_utc = at.Time.now()
                    dt_utc = t_utc.to_datetime(timezone = datetime.timezone.utc)
                    local_timezone = pytz.timezone(timezone_str)
                    dt_local = dt_utc.astimezone(local_timezone)
                    dt_local_str = dt_local.strftime('%Y-%m-%d %H:%M:%S %Z%z')
            
                    return (timezone_str, int(dt_local_str[-5:-2].replace('0', '')))
            
            SN_timezone = [] ; SN_utcoffset = [] ; obs_tz = {}
            tz_offset = {'UTC' : 0, 'Asia/Manila' : 8}
            for lon, lat, elev in zip(SN_long, SN_lat, SN_elev):
                tz, offset = Timezone_Finder(lon, lat, elev)
                SN_timezone.append (tz) ; SN_utcoffset.append (offset)
            for obs, tz, offset in zip(SN_obs, SN_timezone, SN_utcoffset):
                obs_tz.update    ({obs : tz})
                tz_offset.update ({tz  : offset})

            tz_offset = dict(sorted(tz_offset.items(), key=lambda item: item[1]))
            
            obs_list = [SN_obs, SN_long, SN_lat, SN_elev, SN_aperture, SN_telescopes, SN_timezone, SN_utcoffset]
            SN_OBS = pd.DataFrame()
            SN_OBS['Observatory'] = SN_obs
            SN_OBS['Longitude']   = SN_long
            SN_OBS['Latitude']    = SN_lat
            SN_OBS['Elevation']   = SN_elev
            SN_OBS['Timezone']    = SN_timezone
            SN_OBS['UTC Offset']  = SN_utcoffset
            SN_OBS = SN_OBS.sort_values('Observatory')
            st.write ('SkyNet Observatories Geographic Coordinates')
            st.dataframe(SN_OBS)
            # SN_OBS['Aperture']    = SN_aperture
            # SN_OBS['Telescopes']  = SN_telescopes
            # SN_OBS.to_csv('SN_OBS.csv', index = False, header = True)
            
            ####################################################################################################
            st.subheader ('Get_NEAdata')
            st.caption ('This function extracts relevant data of identified exoplanets in the NASA Exoplanet Archive (NEA)')
            
            def Get_NEAdata (targets = None):
                if targets != None:
                    if ',' in targets:
                        targets = targets.split(', ')
                    else:
                        targets = [targets]
                    if type(targets[0]) == list:
                      targets = targets[0]
            
                not_null = []
                select_data = ['hostname'  ,        'pl_name',             'ra',        'dec',
                               'sy_dist'   ,        'sy_snum',        'sy_pnum',
                               'st_mass'   ,         'st_rad',       'pl_masse',    'pl_rade',
                               'pl_eqt'    , #    'pl_eqterr1',     'pl_eqterr2',
                               'st_lum'    ,
                               'pl_orbper' ,  'pl_orbpererr1',  'pl_orbpererr2', 'pl_orbsmax',
                               'pl_tranmid', 'pl_tranmiderr1', 'pl_tranmiderr2',
                               'pl_trandur', 'pl_trandurerr1', 'pl_trandurerr2',
                               'pl_orblper',     'pl_orbtper', 'pl_orbeccen',
                               'pl_trandep'
                              ]
            
                cond_standards = ['hostname' , 'pl_name'   , 'ra', 'dec', 'st_rad', 'pl_rade',
                                  'pl_orbper', 'pl_tranmid', 'pl_trandur']
                not_null    = [column + ' is not null' for column in select_data if column in cond_standards or 'err' in column]
                # not_null    = [column + ' is not null' for column in select_data]
                other_conds = ["discoverymethod = 'Transit'", 'tran_flag = 1']
                where_conds = np.append (not_null, other_conds)
            
                select_text = '' ; where_text = ''
                for d in select_data:
                    if    d == select_data[-1]: select_text += '%s' % d
                    else: select_text += '%s, ' % d
                for c in where_conds:
                    if    c == where_conds[-1]:  where_text += '%s' % c
                    else: where_text += '%s AND ' % c
            
                NEAdata = NEA.query_criteria(
                    table  = "pscomppars",
                    select = select_text,
                    where  = where_text)
                NEAcsv = NEAdata.to_pandas(index = False).sort_values('pl_name')
                if targets != None:
                    NEAcsv = NEAcsv[NEAcsv[['hostname', 'pl_name']].isin(targets).any(axis = 1)]
                NEAcsv = NEAcsv.drop(columns=['sky_coord.ra', 'sky_coord.dec'])
                NEAcsv.rename(columns = {'hostname'  : 'Host Name',
                                          'pl_name'  : 'Planet Name',
                                          'ra'       : 'ra',
                                          'dec'      : 'dec',
                                          'sy_dist'  : 'Distance [pc]',
                                          'sy_snum'  : 'Number of Stars',
                                          'sy_pnum'  : 'Number of Planets',
                                          'st_mass'  : 'Stellar Mass [Solar]',
                                          'st_rad'   : 'Stellar Radius [Solar]',
                                          'st_lum'   : 'Stellar Luminosity [log10(Solar)]',
                                          'pl_masse'        : 'Planet Mass [Earth]',
                                          'pl_rade'         : 'Planet Radius [Earth]',
                                          'pl_eqt'          : 'Planet Temperature [K]',
                                          # 'pl_eqterr1'      : 'Planet Temperature [err 1]',
                                          # 'pl_eqterr2'      : 'Planet Temperature [err 2]',
                                          'pl_trandep'      : 'Transit Depth [%]',
                                          'pl_orbper'       : 'Orbital Period [days]',
                                          'pl_orbpererr1'   : 'Orbital Period [err 1]',
                                          'pl_orbpererr2'   : 'Orbital Period [err 2]',
                                          'pl_orbsmax'      : 'Orbit Semi-Major Axis [au]',
                                          'pl_tranmid'      : 'Transit Midpoint [days]',
                                          'pl_tranmiderr1'  : 'Transit Midpoint [err 1]',
                                          'pl_tranmiderr2'  : 'Transit Midpoint [err 2]',
                                          'pl_trandur'      : 'Transit Duration [hours]',
                                          'pl_trandurerr1'  : 'Transit Duration [err 1]',
                                          'pl_trandurerr2'  : 'Transit Duration [err 2]',
                                          'pl_orblper'      : 'Periastron Argument [deg]',
                                          'pl_orbtper'      : 'Periastron Passage Time [deg]',
                                          'pl_orbeccen'     : 'Eccentricity',
                                          },
                              inplace = True)
                NEAcsv.to_csv('NEAcsv.csv', index = False, header = True)
                NEAcsv = pd.read_csv('NEAcsv.csv')
            
                return (NEAcsv)

            Display_Option = st.radio ('Display the data of all available NEA transiting exoplanets?', ['No', 'Yes'])
            if Display_Option == 'Yes':
                NEAcsv = Get_NEAdata()
                st.dataframe(NEAcsv)
                st.write ('Host Stars:', len(sorted(set(NEAcsv['Host Name']))), '| Exoplanets', len(sorted(set(NEAcsv['Planet Name']))))
            
            ####################################################################################################
            st.subheader('Get_Transits')
            st.caption("This function uses NEA's documented transit algorithm to predict the transit events of target(s), either a specific planet(s) or a system(s), within a specified time window.")
            
            def Get_Transits (targets, start_date, end_date, obs_csv = SN_OBS):
                # targets    = input ('%-64s | ' % "Name of Observation Target:     (Example: Proxima Centauri)")
                if ',' in targets:
                    targets = targets.split(', ')
                elif ',' not in targets and type(targets) != list:
                    targets = [targets]
            
                # start_date = input ('%-64s | ' % "Start Date: YYYY-MM-DDTHH:MM:SS (Example: '2025-12-01T18:00:00')")
                # end_date   = input ('%-64s | ' % "End Date:   YYYY-MM-DDTHH:MM:SS (Example: '2025-12-31T18:00:00')")
                # Day_skip   = input ('%-64s | ' % "Number of day(s) to skip:       (Example: 1)")
            
                # custom_obs = input ('%-64s | ' % "Custom Observatory:             (NAME, LON, LAT, ELEV) / None")
                custom_obs = str(custom_obs)
                if custom_obs != 'None':
                    custom_obs = custom_obs.split(',')
                    lats  = np.append (SN_lat,  float(custom_obs[0])) ; longs = np.append (SN_long, float(custom_obs[1]))
                    elevs = np.append (SN_elev, float(custom_obs[2])) ; obs   = np.append (SN_obs ,       custom_obs[3])
                else:
                    obs = SN_obs ; longs = SN_long ; lats = SN_lat ; elevs = SN_elev
                observatories = []
                observatories.append ([obs, longs, lats, elevs]) ; observatories = observatories[0]
            
                # main_observatory = input ('%-64ss | ' % "Name of Main Observatory:       (Example: Cerro Tololo)")
            
                NEAcsv = Get_NEAdata (targets)
            
                time_tdb = at.Time(start_date).tdb
                time_jd  = time_tdb.jd * u.day
            
                NextTransits_ALL    = []
                NextTransits_Tearly = [] ; NextTransits_err1 = []
                NextTransits_Tlate  = [] ; NextTransits_err2 = []
                end_date  = at.Time(end_date, format = 'iso', scale = 'utc')
                for i in range(len(NEAcsv)):
                    MidTransit           = at.Time(NEAcsv['Transit Midpoint [days]' ][i], format = 'jd', scale = 'utc').jd * u.day
                    MidTransit_err1      =         NEAcsv['Transit Midpoint [err 1]'][i] * u.day
                    MidTransit_err2      =         NEAcsv['Transit Midpoint [err 2]'][i] * u.day
            
                    OrbitalPeriod        = at.Time(NEAcsv['Orbital Period [days]' ]  [i], format = 'jd', scale = 'utc').jd * u.day
                    OrbitalPeriod_err1   =         NEAcsv['Orbital Period [err 1]']  [i] * u.day
                    OrbitalPeriod_err2   =         NEAcsv['Orbital Period [err 2]']  [i] * u.day
            
                    TransitDuration      =         NEAcsv['Transit Duration [hours]'][i] * u.hour
                    TransitDuration_err1 =         NEAcsv['Transit Duration [err 1]'][i] * u.hour
                    TransitDuration_err2 =         NEAcsv['Transit Duration [err 2]'][i] * u.hour
            
                    k = math.ceil((time_jd - MidTransit) / OrbitalPeriod)
                    if k < 0:
                        k = 0
                    epoch_num = 1 ; nexttransits = [0]
                    while float(max(nexttransits)) <= float(end_date.jd):
                        epochs = k + np.arange(epoch_num)
                        nexttransit      = MidTransit      + epochs * np.asarray(OrbitalPeriod     )*u.day
                        NextTransit_err1 = MidTransit_err1 + epochs * np.asarray(OrbitalPeriod_err1)*u.day
                        NextTransit_err2 = MidTransit_err2 + epochs * np.asarray(OrbitalPeriod_err2)*u.day
                        nexttransits = [float(at.Time(jd, format = 'jd', scale = 'utc').value) for jd in nexttransit]
                        epoch_num += 1
            
                    # NextTransits     = [date       for date       in     nexttransits                    if date <= end_date.jd]
                    # # # print (NEAcsv['Planet Name'][i], OrbitalPeriod.value, OrbitalPeriod.value/2)
                    # NextTransits     = [(date + np.round(OrbitalPeriod.value/2)) for date in nexttransits          if date <= end_date.jd]
                    NextTransits     = [date for date in nexttransits          if date <= end_date.jd]
                    NextTransit_err1 = [err1.value for date, err1 in zip(nexttransits, NextTransit_err1) if date <= end_date.jd]
                    NextTransit_err2 = [err2.value for date, err2 in zip(nexttransits, NextTransit_err2) if date <= end_date.jd]
            
                    T_early = epochs * np.asarray(OrbitalPeriod_err1)*u.day + MidTransit_err1 + 0.5 * TransitDuration + TransitDuration_err1
                    T_late  = epochs * np.asarray(OrbitalPeriod_err2)*u.day + MidTransit_err2 + 0.5 * TransitDuration + TransitDuration_err2
            
                    NextTransits_ALL.append    (NextTransits)
                    NextTransits_Tearly.append ((T_early/u.day).value) ; NextTransits_err1.append (NextTransit_err1)
                    NextTransits_Tlate.append  ((T_late /u.day).value) ; NextTransits_err2.append (NextTransit_err2)
            
                NEAcsv['Next Transits [JD]'   ] = NextTransits_ALL
                NEAcsv['Next Transits [err 1]'] = NextTransits_err1   ; NEAcsv['Next Transits [err 2]'] = NextTransits_err2
                NEAcsv['Next Transits [early]'] = NextTransits_Tearly ; NEAcsv['Next Transits [late]' ] = NextTransits_Tlate
            
                NEAcsv.to_csv('NEAcsv.csv', index = False, header = True)
            
                return (NEAcsv)

            NEAcsv = Get_Transits(
            targets    = targets,
            start_date = str(start_date),
            end_date   = str(end_date),
            )
            
            # NEAcsv = NEAcsv[NEAcsv['Host Name'].isin(['TRAPPIST-1', 'Kepler', 'WASP', 'XO'])].reset_index()
            st.dataframe(NEAcsv)
            
            ####################################################################################################
            st.subheader ('Generate_Transit_Dates')
            st.caption   ("This function extracts only the visible transit events per observatory in SkyNet.")
            st.caption   ("The data in this table are the ingress, midpoint, and the egress of the predicted transit events respectively.")
            
            ####################################################################################################
            st.subheader ('Generate_Transit_Dates')
            st.caption   ("This function extracts only the visible transit events per observatory in SkyNet.")
            st.caption   ("The data in this table are the ingress, midpoint, and the egress of the predicted transit events respectively.")
            
            def Generate_Transit_Dates (CSV, timezone = 'UTC', obs_csv = SN_OBS):
                st.write('Generate_Transit_Dates || %s %s' % (timezone, tz_offset[timezone]))
                obs_names = obs_csv['Observatory'] ; longs = obs_csv['Longitude'] ; lats = obs_csv['Latitude'] ; elevs = obs_csv['Elevation']
                DF = pd.DataFrame()
                Obs_All     = [] ; Host_All     = [] ; Planet_All = []
                Ingress_All = [] ; Midpoint_All = [] ; Egress_All = []
            
                constraints = [ap.AltitudeConstraint(min = 20*u.deg, max = 90*u.deg),
                               ap.AtNightConstraint(-12*u.deg)]
            
                for idx in range(len(CSV)):
                    P_Ingress = [] ; P_Midpoint = [] ; P_Egress = []
                    ra = CSV['ra'][idx] ; dec = CSV['dec'][idx]
                    HostName = CSV['Host Name'][idx] ; PlanetName = CSV['Planet Name'][idx]
                    NextTransits_ALL = NEAcsv['Next Transits [JD]'   ][idx]
                    Tearly_ALL       = NEAcsv['Next Transits [early]'][idx] ; Tlate_ALL = NEAcsv['Next Transits [late]'][idx]
                    for obs_name, lon, lat, elev in zip(obs_names, longs, lats, elevs):
                        Obs_Ingress = [] ; Obs_Midpoint = [] ; Obs_Egress = []
                        Host_All.append (HostName) ; Planet_All.append (PlanetName) ; Obs_All.append (obs_name)
                        location = ac.EarthLocation.from_geodetic (lon, lat, elev*u.m)
                        observer = ap.Observer (location = location, name = obs_name)
                        coord    = SkyCoord (ra = ra, dec = dec, unit = 'deg')
                        target   = FixedTarget (coord = coord, name = obs_name)
            
                        for Tmid, Tearly, Tlate in zip(NextTransits_ALL, Tearly_ALL, Tlate_ALL):
                            T_mid      = at.Time(Tmid         , format = 'jd', scale = 'utc')
                            Day_Tearly = at.Time(Tmid - Tearly, format = 'jd', scale = 'utc')
                            Day_Tlate  = at.Time(Tmid - Tlate , format = 'jd', scale = 'utc')
                            if ap.is_observable(constraints, observer, target, time_range = [Day_Tearly, Day_Tlate]) == True:
                                T_mid      += tz_offset[timezone]*u.hour
                                Day_Tearly += tz_offset[timezone]*u.hour
                                Day_Tlate  += tz_offset[timezone]*u.hour
                                Obs_Ingress.append (Day_Tearly.iso) ; Obs_Midpoint.append (T_mid.iso) ; Obs_Egress.append (Day_Tlate.iso)
                        Ingress_All.append (Obs_Ingress) ; Midpoint_All.append (Obs_Midpoint) ; Egress_All.append (Obs_Egress)
            
                DF['Host Name'  ] = Host_All
                DF['Planet Name'] = Planet_All
                DF['Observatory'] = Obs_All
                DF['Ingress'    ] = Ingress_All
                DF['Midpoint'   ] = Midpoint_All
                DF['Egress'     ] = Egress_All
            
                DF.to_csv ('Observatory Visible Transit Dates - UTC.csv', index = False, header = True)
                DF_utc = DF
            
                return (DF_utc)
            
            options  = [(tz, offset) for tz, offset in tz_offset.items()]
            timezone = st.selectbox (label = 'Timezone of output dates / UTC offset', options = options, index = options.index(('UTC', 0)))
            tz_cb    = st.checkbox  ('Submit Timezone')
            if tz_cb:
                TDates = Generate_Transit_Dates (NEAcsv, timezone, SN_OBS)
                st.dataframe(TDates)
                
                ####################################################################################################
                # st.subheader('Generate_Calendars')
                # st.caption("This function produces an .ics file for users to import in their calendar (Gmail and/or Outlook Calendar) for their planning purposes.")
                
                # # @st.cache_data
                # def Generate_Calendar (TDatesCSV):
                #     c = Calendar()
                #     obs_name = TDates['Observatory']
                #     for PlanetName in TDates:
                #         if PlanetName == 'Observatory': continue
                #         for idx, obs_name in enumerate (TDates['Observatory']):
                #             for Dates in TDates[PlanetName]:
                #                 for Date in Dates:
                #                     if len(Dates) == 0: continue
                #                     e = Event()
                #                     e.name = PlanetName
                #                     e.begin = at.Time(Date[0]).datetime
                #                     e.end   = at.Time(Date[1]).datetime
                #                     e.location = obs_name
                #                     c.events.add(e)
                #     with open('TransitDates.ics', 'w') as f:
                #         f.writelines(c.serialize_iter())
                #     st.download_button(label = "Download ICS file",
                #                        data = f,
                #                        file_name = "TransitDates.ics",
                #                        mime="text/csv",
                #                        icon=":material/download:",
                #                       )
                # Generate_Calendar (TDates)
    
                ####################################################################################################
                # def Generate_Transit_Airmass (CSV, obs_list, custom_obs = None, main_observatory = 'Cerro Tololo'):
                #     if os.path.isdir(os.getcwd() + '/Transits') == False:
                #         os.mkdir (os.getcwd() + '/Transits')
                
                #     Err_list = []
                #     mpl.rcParams['axes.linewidth'] = 4
                #     min_alt = 20 ; max_alt = 90
                #     label_fs = 16
                #     lwidth = 2.25
                
                #     hour = 18
                #     hour_offset  = 1 - hour/24
                #     constraints = [ap.AltitudeConstraint(min = 20*u.deg, max = 90*u.deg),
                #                    ap.AtNightConstraint(-12*u.deg)]
                
                #     for idx in range(len(CSV)):
                #         starttime = timeit.default_timer()
                #         host = CSV['Host Name'][idx] ; planet = CSV['Planet Name'][idx]
                #         if os.path.isdir (os.getcwd() + '/Transits/%s' % host) == False:
                #             os.mkdir     (os.getcwd() + '/Transits/%s' % host)
                
                #         ra   = CSV['ra'][idx] ; dec = CSV['dec'][idx] ; Err_list = []
                #         Days = CSV['Next Transits [JD]'][idx]
                
                #         observables = 0
                #         for day_count, Day in enumerate(Days):
                #             try:
                #                 night_check = 1
                #                 Day_iso = at.Time(Day, format = 'jd', scale = 'utc').iso
                #                 td_start = Day
                #                 if np.floor(Day) + hour_offset < Day and Day <= np.ceil(Day) + hour_offset:
                #                     td_start = np.floor(Day) + hour_offset
                #                 elif np.floor(Day) + hour_offset > Day:
                #                     td_start = np.floor(Day) - (1 - hour_offset)
                #                 td_start     = td_start * u.day + np.linspace (0, +24, 97) * u.hour
                #                 time         =     at.Time(td_start, format = 'jd', scale = 'utc').iso
                #                 current_date = str(at.Time(     Day, format = 'jd', scale = 'utc').iso).split(' ')[0]
                
                #                 NextTransit_Tearly = CSV['Next Transits [early]'][idx][day_count] ; NextTransit_err1 = CSV['Next Transits [err 1]'][idx][day_count]
                #                 NextTransit_Tlate  = CSV['Next Transits [late]' ][idx][day_count] ; NextTransit_err2 = CSV['Next Transits [err 2]'][idx][day_count]
                #                 Day_Tearly = at.Time (Day - NextTransit_Tearly              , format = 'jd', scale = 'utc')
                #                 Day_Tlate  = at.Time (Day - NextTransit_Tlate               , format = 'jd', scale = 'utc')
                #                 Day_T_mid  = at.Time ((Day_Tearly.value + Day_Tlate.value)/2, format = 'jd', scale = 'utc')
                #                 # Day_T_mid  = at.Time (Day                                   , format = 'jd', scale = 'utc')
                #                 Day_Terr1  = at.Time (Day - NextTransit_err1                , format = 'jd', scale = 'utc')
                #                 Day_Terr2  = at.Time (Day - NextTransit_err2                , format = 'jd', scale = 'utc')
                
                #                 for obs_name, lon, lat, elev in zip(obs_list[0], obs_list[1], obs_list[2], obs_list[3]):
                #                     if main_observatory in obs_name:
                #                         location = ac.EarthLocation.from_geodetic (lon, lat, elev*u.m)
                #                         observer = ap.Observer (location = location, name = obs_name)
                #                         coord    = SkyCoord (ra = ra, dec = dec, unit = 'deg')
                #                         target   = FixedTarget (coord = coord, name = obs_name)
                #                         if ap.is_observable(constraints, observer, target, time_range=[Day_Tearly, Day_Tlate]) == False:
                #                             night_check = 0
                #                             break
                #                 if night_check == 0: continue
                
                #                 if os.path.isdir (os.getcwd() + '/Transits/%s/%s' % (host, planet)) == False:
                #                     os.mkdir     (os.getcwd() + '/Transits/%s/%s' % (host, planet))
                
                #                 observables += 1
                #                 plt.figure(figsize = (24, 12))
                #                 for obs_name, lon, lat, elev in zip(obs_list[0], obs_list[1], obs_list[2], obs_list[3]):
                #                     location = ac.EarthLocation.from_geodetic (lon, lat, elev*u.m)
                #                     observer = ap.Observer (location = location, name = obs_name)
                
                #                     coord  = SkyCoord (ra = ra, dec = dec, unit = 'deg')
                #                     target = FixedTarget (coord = coord, name = obs_name)
                
                #                     # if observer.target_is_up(at.Time(Day, format = 'jd', scale = 'utc').iso, target, horizon = 20*u.deg) == False:
                #                     if ap.is_observable(constraints, observer, target, time_range=[Day_Tearly, Day_Tlate]) == False:
                #                         continue
                
                #                     # plotting
                #                     msize = 0
                #                     ls = {'linewidth' : lwidth, 'markersize' : msize}
                #                     if lat <  0: ls.update({'linestyle' : '--'})
                #                     if lat >= 0: ls.update({'linestyle' : '-' })
                #                     # if 'Green Bank' in obs_name:
                #                     #     ls.update({'linestyle' : '-.'})
                #                     if main_observatory in obs_name:
                #                         m_ls    = {'linewidth' : lwidth * 2.35, 'linestyle' : ':' , 'markersize' :  msize, 'color' : 'black'}
                #                         ls.update({'linewidth' : lwidth * 2.35, 'color' : 'black'})
                #                         moon = observer.moon_altaz (time) ; moon_radec = moon.transform_to('icrs')
                #                         moon_ra = moon_radec.ra ; moon_dec = moon_radec.dec
                #                         moon_coord  = SkyCoord (ra = moon_ra, dec = moon_dec, unit ='deg', frame = 'icrs')
                #                         moon_target = FixedTarget (coord = moon_coord, name = 'Moon (%s)' % main_observatory)
                #                         ap.plots.plot_altitude (     target, observer, time,
                #                                                 brightness_shading = True,  min_altitude = min_alt, max_altitude = max_alt, style_kwargs =   ls)
                #                         ap.plots.plot_altitude (moon_target, observer, time,
                #                                                 brightness_shading = False, min_altitude = min_alt, max_altitude = max_alt, style_kwargs = m_ls)
                #                         AltAz_frame1  = AltAz (obstime = Day_Tearly, location = location)
                #                         AltAz_frame2  = AltAz (obstime = Day_Tlate , location = location)
                #                         AltAz_frameM  = AltAz (obstime = Day_T_mid , location = location)
                #                         target_AltAz1 = coord.transform_to(AltAz_frame1) ; ingress_alt = target_AltAz1.alt.value
                #                         target_AltAz2 = coord.transform_to(AltAz_frame2) ; egress_alt  = target_AltAz2.alt.value
                #                         target_AltAzM = coord.transform_to(AltAz_frameM) ; mid_alt     = target_AltAzM.alt.value
                
                #                     else:
                #                         ap.plots.plot_altitude (     target, observer, time,
                #                                                 brightness_shading = False, min_altitude = min_alt, max_altitude = max_alt, style_kwargs =   ls)
                
                #                 plt.xticks(np.linspace(plt.xticks()[0][0], plt.xticks()[0][-1], 13).tolist())
                #                 # xlabel_date = plt.gca().get_xlabel().split(' ')[2]
                #                 plt.title ('%s on %s' % (planet, current_date), fontsize = label_fs)
                #                 raw_xticks = [x for x in plt.gca().get_xticks()]
                #                 MIN = min(raw_xticks) ; MAX = max(raw_xticks)
                #                 Day_iso = at.Time(Day, format = 'jd', scale = 'utc').iso
                
                #                 plt.axvspan (Day_Tearly.datetime, Day_Tlate.datetime, color = 'green', alpha = 0.3) # ingress, egress
                #                 plt.axvline ( Day_T_mid.datetime, linewidth = 1.5   , color = 'green', alpha = 0.3) # window midpoint
                #                 # print ('\n', Day_Terr1.datetime, '|', Day_Terr2.datetime, '|', (Day_Terr2 - Day_Terr1).value*u.day.to(u.hour), (Day_Terr2 - Day_Terr1).value*u.day.to(u.minute))
                #                 # plt.axvline ( Day_Terr1.datetime, linewidth = 1.5   , color = 'orange', alpha = 0.3)
                #                 # plt.axvline ( Day_Terr2.datetime, linewidth = 1.5   , color = 'orange', alpha = 0.3)
                #                 plt.text (Day_T_mid.datetime, mid_alt, fontsize = label_fs,
                #                           s = '%+8s – %-8s\n %+10s – %s' % (str(Day_Tearly.datetime)[10:19], str(Day_Tlate.datetime)[10:19],
                #                                                             str(np.round(ingress_alt, 2)) + '\N{DEGREE SIGN}',
                #                                                             str(np.round(egress_alt , 2)) + '\N{DEGREE SIGN}'))
                
                #                 ######## LABELS
                #                 #### x-labels
                #                 plt.xlabel('Time [UTC]', size = label_fs)
                #                 # plt.xlabel('', size = label_fs + 4)
                #                 zero_UTC  = [i for i, x in enumerate(raw_xticks) if x == int(x)][0]
                #                 plt.gca().get_xticklabels()[zero_UTC].set_weight ("bold")
                #                 plt.gca().get_xticklabels()[zero_UTC].set_color  ("red")
                #                 #### y-labels
                #                 plt.ylabel('Altitude $[\u00b0]$', size = label_fs)
                #                 # plt.ylabel('', size = label_fs + 4)
                #                 ay2 = plt.gca().secondary_yaxis('right')
                #                 custom_locs   = np.linspace(10, 90, 9, dtype = int).tolist()
                #                 custom_yticks = [5.86, 2.9, 1.993, 1.553, 1.304, 1.154, 1.064, 1.015, 1.0]
                #                 ay2.set_yticks (custom_locs)
                #                 ay2.set_yticklabels (custom_yticks, fontsize = label_fs)
                #                 ay2.set_ylabel ('Airmass', rotation = 270, labelpad = 18, size = label_fs)
                #                 # ay2.set_ylabel ('', rotation = 270, labelpad = 18, size = label_fs)
                
                #                 plt.tick_params (axis = 'both', labelsize = label_fs)
                #                 plt.xticks (fontsize = label_fs * 1.5)
                #                 plt.yticks (fontsize = label_fs * 2)
                
                #                 plt.grid (axis = 'y', linewidth = 3.25, alpha = 0.5)
                
                #                 first_legend = plt.legend (fontsize = label_fs, bbox_to_anchor = (1, 1),
                #                                             loc = 'upper right', title = 'Visible', title_fontsize = label_fs, framealpha = 0.9)
                #                 custom_line1 = mlines.Line2D ([], [], c = 'dimgrey', ls =  '-', lw = lwidth, label = 'Northern Hemisphere')
                #                 custom_line2 = mlines.Line2D ([], [], c = 'dimgrey', ls = '--', lw = lwidth, label = 'Southern Hemisphere')
                #                 secondary_legend = plt.legend (fontsize = label_fs, handles=[custom_line1, custom_line2], loc='upper left', frameon = False)
                #                 plt.gca().add_artist(secondary_legend)
                #                 plt.gca().add_artist(first_legend)
                
                #                 plt.gcf().set_size_inches (24, 12)
                #                 filename = os.getcwd() +'/Transits/%s/%s/AIR_%s.png' % (host, planet, current_date)
                #                 plt.savefig (filename, bbox_inches = 'tight', dpi = 600)
                #                 plt.close ('all')
                #                 print ('\r Progress: %+15s | %+15s | %-10s | %+3s/%-3s days' %
                #                     (host, planet, current_date, day_count, len(Days)), end = '')
                #             except Exception as err:
                #                 Err_list.append (err)
                
                #         stoptime = timeit.default_timer()
                #         runtime  = stoptime-starttime ; minute = math.floor(runtime/60) ; seconds = runtime - 60*minute
                #         if observables == 0:   print (  '      | %s has no visible transit.' % planet)
                #         else:                  print ("\rDONE! | Runtime: %+4sm %+8s %-50s" % (int(minute), str(np.round(seconds, 2)) + 's |', (host + ' — ' + planet)))
                #         if len(Err_list) != 0: print ('Errors:\n', Err_list)
