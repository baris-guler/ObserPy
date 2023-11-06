from datetime import datetime, timedelta, date, time, timezone
import datetime as dt
import math
import numpy as np
import pandas
import re
from PyQt5 import QtWidgets, uic
from PyQt5.QtWidgets import QDialog, QVBoxLayout, QLabel, QLineEdit, QPushButton, QTableWidgetItem, QFileDialog, QMessageBox, QHeaderView, QCheckBox, QTextEdit
from decimal import *
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as NavigationToolbar)
from PyQt5.QtGui import QColor, QStandardItem, QStandardItemModel
import ephem
from PyQt5 import QtCore
from PyQt5.QtCore import Qt, QRunnable, QThreadPool, pyqtSignal, QObject, pyqtSlot, QDateTime, QDate, QTime
import matplotlib.dates as mdates
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery import simbad
import matplotlib.lines as mlines
from astroquery.mast import Observations
from astropy.io import fits
import os
from astroquery.mast import Tesscut
import pandas as pd
import webbrowser
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
import pickle
import log
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import csv

class Location:
    def __init__(self, name, long, lat_deg, tbz_gmt):
        self.name = name
        self.long = long
        self.long_h = float(long)/15
        self.lat_deg = float(lat_deg)
        self.tbz_gmt = int(tbz_gmt)

class Star:
    def __init__(self, name, ra_h, dec_deg, E, P, mag, pmra=None, pmdec=None, min_start=np.nan, min_end=None):
        for i in stars:
            if i.name == name:
                raise IndexError("This Star Already Added:", name)
        self.name = name
        if (pmra is None) | (pmdec is None):
            self.pmra, self.pmdec = 0,0
        else:
            self.pmra, self.pmdec = pmra, pmdec
        self.ra = float(self._hours_to_h(ra_h))
        self.dec = float(self._deg_to_float(dec_deg))
        self.E = float(E)
        self.P = float(P)
        self.h = []
        self.mag = float(mag)
        self.mintimes = []
        self.mintimes2 = []
        self.minstart = min_start
     
    @staticmethod
    def get_proper_motion_from_simbad(object_name):
        custom_simbad = simbad.Simbad()
        custom_simbad.add_votable_fields('pmra', 'pmdec')
        result_table = custom_simbad.query_object(object_name)

        if result_table is not None and len(result_table) > 0:
            pmra = result_table['PMRA'][0]
            pmdec = result_table['PMDEC'][0]
            return pmra / 3600000.0, pmdec / 3600000.0
        else:
            return 0, 0
    
    def _hours_to_h(self, ra_h):
        global start_date
        try:
            ra_h = float(ra_h)
        except:
            pass
        if isinstance(ra_h, str):
            self.ra_str = ra_h
            ra_h = ra_h.strip()
            if (ra_h.find(":") != -1):
                splitter = ":"
            else:
                ra_h = re.sub('\s+', ' ', ra_h)
                splitter = " "
            h, min, sec = ra_h.split(splitter)
            hours_float = (int(h) + int(min)/60 + float(sec)/3600)
            return hours_float/24 + (start_date.year-2000)*self.pmra/360
        else:
            return ra_h + (start_date.year-2000)*self.pmra/360

    def _deg_to_float(self, lat_deg):
        global start_date
        try:
            lat_deg = float(lat_deg)
        except:
            pass
        if isinstance(lat_deg, str):
            self.dec_str = lat_deg
            lat_deg = lat_deg.strip()
            if (lat_deg.find(":") != -1):
                splitter = ":"
            else:
                lat_deg = re.sub('\s+', ' ', lat_deg)
                splitter = " "
            deg, arcmin, arcsec = lat_deg.split(splitter)
            if deg[0] == "-":
                sign = -1
            else:
                sign = 1
            hours_float = (int(deg) + int(arcmin)/60*sign + float(arcsec)/3600*sign)
            return hours_float + (start_date.year-2000)*self.pmra
        else:
            return lat_deg + (start_date.year-2000)*self.pmra
        
class WorkerSignals(QObject):
    finished = pyqtSignal(object)  # Add an object signal

class Worker(QRunnable):
    def __init__(self, function):
        super().__init__()
        self.function = function
        self.signals = WorkerSignals()

    @pyqtSlot()
    def run(self):
        try:
            result = self.function()  # Get the return value of the function
            self.signals.finished.emit(result)  # Emit the result
        except:
            pass
        
class AboutDialog(QDialog):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("About My App")
        self.setGeometry(500, 300, 800, 300)

        layout = QVBoxLayout()

        about_text = """
        ---------------------------------------------------------------------------------------------------------------------------------------------
        Selecting a specific system for the observations of extrema in the light curves of variable stars poses a significant challenge, 
        particularly in cases where observers are confronted with a multitude of systems but face time constraints. In response to this issue, 
        we developed an application with Python to help the observer to determine the most viable system to be observed at any given time. 
        This application also provides information on optimal dates for observing a particular system at a specified time interval. Additionally, 
        the application enables users to create detailed logs documenting any observation of a particular target.
        This is a free, open source program.
        ---------------------------------------------------------------------------------------------------------------------------------------------
        Verison: 1.0
        Github Link: https: github.com/baris-guler/ObservaPy
        Author: Barış Güler
        My mail: barisguler2000@gmail.com
        """

        text_edit = QTextEdit()
        text_edit.setPlainText(about_text)
        text_edit.setReadOnly(True)  # Make the text uneditable
        layout.addWidget(text_edit)

        self.setLayout(layout)
        
def get_julian_datetime(date):
    # Ensure correct format
    if not isinstance(date, dt.datetime):
        raise TypeError('Invalid type for parameter "date" - expecting datetime')
    elif date.year < 1801 or date.year > 2099:
        raise ValueError('Datetime must be between year 1801 and 2099')

    # Perform the calculation
    date = datetime(date.year, date.month, date.day)
    julian_datetime = 367 * date.year - int((7 * (date.year + int((date.month + 9) / 12.0))) / 4.0) + int((275 * date.month) / 9.0) + date.day + 1721013.5 + (date.hour + date.minute / 60.0 + date.second / math.pow(60, 2)) / 24.0 - 0.5 * math.copysign(1, 100 * date.year + date.month - 190002.5) + 0.5

    return julian_datetime - 2400000

def hours_to_days(hour=0, minute=0, seconds=0):
    d = hour + minute/60 + seconds/3600
    return d/24

def load():
    global obs_loc
    global star_file
    global observatory_file
    global stars
    global df
    global observatories
    
    dlg.star_table.clearContents()
    df = pd.DataFrame()
    stars = []
    observatories = []
    if star_file.endswith(".csv"):
        with open(star_file, "r") as file:
            first_line = file.readline()
            if "," in first_line:
                loaded_data = pd.read_csv(star_file)
            else:
                loaded_data = pd.read_csv(star_file, delimiter=";")
        for index, row in loaded_data.iterrows():
            sname = row["Star Name"]
            sra = row["ra(h)"]
            sdec = row["dec(deg)"]
            sE = row["Epoch"]
            sP = row["Period"]
            smag = row["Magnitude"]
            mintime_start = row["Mintime_Start"]
            stars.append(Star(sname, sra, sdec, sE, sP, smag, min_start=mintime_start))
    else:
        raise ValueError("Star file must be .txt or .csv")
            
    if observatory_file.endswith(".csv"):
        obs_df = pd.read_csv(observatory_file, sep=",")
        for index, row in obs_df.iterrows():
            oname = row["Observatory"]
            olong = row["Longitude"]
            olat = row["Latitude"]
            otbz_gmt = row["TBZ_GMT"]
            observatories.append(Location(oname, olong, olat, otbz_gmt))
        obs_loc = observatories[0]
    else:
        raise ValueError("Observatory file must be .txt or .csv")


def gst_jd_calculator(time):
    if time.hour < 12:
        jd_ut = get_julian_datetime(time) - 1
    else:
        jd_ut = get_julian_datetime(time)
    T = (jd_ut-15020)/36525
    gst_ut = (((0.276919398)+100.0021359*T+0.000001075*T**2)-int(0.276919398+100.0021359*T+0.000001075*T**2))*24
    
    tbz = time
    ut = tbz - timedelta(hours = obs_loc.tbz_gmt)
    if ut.hour >= 0 and ut.hour < 9:
        mult = 1
    else:
        mult = 0
    gst = gst_ut + (hours_to_days(ut.hour+24*mult, ut.minute, ut.second)*24*1.002737908)
    jd = jd_ut + hours_to_days(ut.hour+24*mult + ut.minute/60 + ut.second/3600)
    return gst, jd

def h_calculator(t, star):
    global s
    global times
    global min_alt
    global M
    
    gst, jd = gst_jd_calculator(t)
    LST = (gst - obs_loc.long_h)/24
    
    if math.copysign(1, LST - star.ra) == -1:
        HA = abs((LST-star.ra) - math.trunc(LST-star.ra))
    elif ((LST-star.ra)-math.trunc(LST-star.ra))*24>12:
        HA = 1-((LST-star.ra) - math.trunc(LST-star.ra))
    else:
        HA = (LST-star.ra) - math.trunc(LST-star.ra)
    if math.copysign(1, LST - star.ra) == -1:
        sign_HA = -1
    elif ((LST - star.ra)-math.trunc(LST - star.ra))*24>12:
        sign_HA = -1
    else:
        sign_HA = 1
    HA = abs((LST-star.ra) - math.trunc(LST-star.ra))
    h = math.degrees(math.asin(math.sin(math.radians(obs_loc.lat_deg))*math.sin(math.radians(star.dec))+math.cos(math.radians(obs_loc.lat_deg))*
                            math.cos(math.radians(star.dec))*math.cos(math.radians(HA*15*sign_HA*24))))
    return h
            
def dataframe_refresh():
    global s
    global min_alt
    global M
    global twilight_start
    global twilight_end
    global observer
    global df
    global columns_active
    if "obs" in df.columns:
        df.drop("obs", axis=1, inplace=True)
    if "_mintime_start_h" in df.columns:
        df.drop("_mintime_start_h", axis=1, inplace=True)
    
    #Moon refrsesh
    observer.date = f"{start_date.year}/{start_date.month}/{start_date.day} {start_date.hour-obs_loc.tbz_gmt}:{start_date.minute}" 
    M = ephem.Moon(observer)
    M.compute(observer)
    stars[-1].ra = stars[-1]._hours_to_h(str(M.ra))
    stars[-1].dec = stars[-1]._deg_to_float(str(M.dec))
    moon_skycoord = SkyCoord(ra=stars[-1].ra*360*u.deg, dec=stars[-1].dec*u.deg, frame='icrs')
    
    fred = ephem.Observer()

    if start_date.hour < 12:
        fred.date = start_date.replace(hour=12) - timedelta(days=1)
    else:
        fred.date = start_date.replace(hour=12)

    fred.lon  = str(obs_loc.long_h*15*-1) 
    fred.lat  = str(obs_loc.lat_deg)      

    fred.elev = 20

    fred.pressure= 0
    fred.horizon = '-18'
    twilight_end=fred.next_rising(ephem.Sun(), use_center=True).datetime() + timedelta(hours=obs_loc.tbz_gmt)
    twilight_start=fred.next_setting(ephem.Sun(), use_center=True).datetime() + timedelta(hours=obs_loc.tbz_gmt)
    
    times = []
    if start_date > twilight_start:
        current_time = start_date
    else:
        current_time = twilight_start
        
    while current_time < twilight_end:
        times.append(current_time)
        current_time += timedelta(minutes=1)
    
    for star in stars:
        gst, jd = gst_jd_calculator(start_date)
        h = h_calculator(start_date, star)
        
        #mintimes
        star.mintimes = []
        star.countable_mintimes = []
        time = start_date
        ph = ((jd-star.E)/star.P)%1
        while time > twilight_start - timedelta(hours=5):
            time -= timedelta(hours=(ph)*star.P*24)
            ph = 1
            if time > twilight_start  - timedelta(hours=5):
                star.mintimes.append(time)
        ph = ((jd-star.E)/star.P)%1
        time = start_date 
        while time < twilight_end + timedelta(hours=5):
            time += timedelta(hours=(1-ph)*star.P*24)
            ph = 0
            if time < twilight_end + timedelta(hours=5) and time > twilight_start - timedelta(hours=5):
                star.mintimes.append(time)
                    
        #mintimes2
        star.mintimes2 = []
        star.countable_mintimes2 = []
        time = start_date
        ph = ((jd-star.E)/star.P)%1 + .5
        while time > twilight_start - timedelta(hours=5):
            time -= timedelta(hours=(ph)*star.P*24)
            ph = 1
            if time > twilight_start - timedelta(hours=5):
                star.mintimes2.append(time)
        ph = ((jd-star.E)/star.P)%1 + .5
        time = start_date  
        while time < twilight_end + timedelta(hours=5):
            time += timedelta(hours=(1-ph)*star.P*24)
            ph = 0
            if time < twilight_end + timedelta(hours=5) and time > twilight_start - timedelta(hours=5):
                star.mintimes2.append(time)
        
        
        if mst != 0:
            for m in star.mintimes:
                if np.isnan(star.minstart):
                    mintime_9 = m - timedelta(minutes=mst)
                else:
                    mintime_9 = m - timedelta(minutes=star.minstart)
                if mintime_9 > start_date:
                    if h_calculator(mintime_9, star) > min_alt:
                        mintime = m
                        if h_calculator(mintime, star) > min_alt:
                            if np.isnan(star.minstart):
                                mintime_11 = m + timedelta(minutes=mst)
                            else:
                                mintime_11 = m + timedelta(minutes=star.minstart)
                            if mintime_11 < twilight_end and h_calculator(mintime_11, star) > min_alt and mintime_9 > twilight_start:
                                observable = True
                                break
                                
            else:
                for m in star.mintimes:
                    if mintime_9 > start_date:
                        if np.isnan(star.minstart):
                            mintime_9 = m - timedelta(minutes=mst)
                        else:
                            mintime_9 = m - timedelta(minutes=star.minstart)
                        mintime = m
                        if np.isnan(star.minstart):
                            mintime_11 = m + timedelta(minutes=mst)
                        else:
                            mintime_11 = m + timedelta(minutes=star.minstart)
                        observable = False
                        break
                else:
                    mintime_9 = start_date
                    mintime = start_date
                    mintime_11 = start_date
                    observable = False
        else:
            for m in star.mintimes:
                mintime_9 = m - timedelta(minutes=mst)
                mintime = m
                mintime_11 = m + timedelta(minutes=mst)
                
            for t in times:
                if h_calculator(t, star) > min_alt:
                    observable = True
                    break
            else:
                observable = False
            
                
        ph = ((jd-star.E)/star.P)%1
        #h = h_calculator(start_date, star)
                
        
        if star.name == "moon":
            dlg.moon_phase.setText(f"%{M.phase:.2f}")    
        else:
            #dataframe
            
            df.loc[star.name, "Star Name"] = star.name
            if columns_active[0]:
                df.loc[star.name, "Cur Altitude"] = round(h, 3)
            elif "Cur Altitude" in df.columns:
                df.drop("Cur Altitude", axis=1, inplace=True)
            if columns_active[1]:
                df.loc[star.name, "Cur Phase"] = round(ph, 3)
            elif "Cur Phase" in df.columns:
                df.drop("Cur Phase", axis=1, inplace=True)
            min_time_start = mintime_9
            min_time_min = mintime
            min_time_end = mintime_11
            if columns_active[2]:
                df.loc[star.name, "Period"] = f"{star.P:.4f}"
            elif "Period" in df.columns:
                df.drop("Period", axis=1, inplace=True)
            if columns_active[3]:
                df.loc[star.name, "Start"] = f"{min_time_start.hour:02d}:{min_time_start.minute:02d}"
            elif "Start" in df.columns:
                df.drop("Start", axis=1, inplace=True)
            if columns_active[4]:
                df.loc[star.name, "End"] = f"{min_time_end.hour:02d}:{min_time_end.minute:02d}"
            elif "End" in df.columns:
                df.drop("End", axis=1, inplace=True)
            if columns_active[5]:
                df.loc[star.name, "Mintime"] = f"{min_time_min.hour:02d}:{min_time_min.minute:02d}"
            elif "Mintime" in df.columns:
                df.drop("Mintime", axis=1, inplace=True)
            if columns_active[6]:
                df.loc[star.name, "Mag"] = star.mag
            elif "Mag" in df.columns:
                df.drop("Mag", axis=1, inplace=True)
            
            star_skycoord = SkyCoord(ra=star.ra*360*u.deg, dec=star.dec*u.deg, frame='icrs')
            if columns_active[7]:
                df.loc[star.name, "Mintime Alt"] = int(round(h_calculator(min_time_min, star)))
            elif "Mintime Alt" in df.columns:
                df.drop("Mintime Alt", axis=1, inplace=True)
            if columns_active[8]:
                df.loc[star.name, "Moon Distance"] = int(str(star_skycoord.separation(moon_skycoord)).split("d")[0])
            elif "Moon Distance" in df.columns:
                df.drop("Moon Distance", axis=1, inplace=True)
            
            if observable:
                df.loc[star.name, "obs"] = True
            else:
                df.loc[star.name, "obs"] = False
            df.loc[star.name, "_mintime_start_h"] = round((min_time_start.day - start_date.day + (min_time_start.hour - start_date.hour)/24 + (min_time_start.minute - start_date.minute)/24/60), 4)
          
    
#GUI
def add_star_dialog():
    dialog = QDialog()
    dialog_layout = QVBoxLayout()
    dialog.setWindowTitle('Add New Star')
    dialog.closeEvent = close_dialog
    label1 = QLabel('Name')
    text1 = QLineEdit("GK Vir")
    label2 = QLabel('RA')
    text2 = QLineEdit("+50 13 10.62")
    label3 = QLabel('Dec')
    text3 = QLineEdit("1:17:18")
    label4 = QLabel('E')
    text4 = QLineEdit("41220.315")
    label5 = QLabel('P')
    text5 = QLineEdit("0.23369101")
    label6 = QLabel('M')
    text6 = QLineEdit("10.62")
    label8 = QLabel('Mintime - Mintime Start (minute)(optional)')
    text8 = QLineEdit("")
    label7 = QLabel()
    add_button = QPushButton("add")
    add_button.clicked.connect(lambda: new_star_add(text1.text(), text2.text(), text3.text(), text4.text(), text5.text(), text6.text(), text8.text(), label7, dialog))
    back_button = QPushButton("back")
    back_button.clicked.connect(lambda: back(dialog))
    
    dialog_layout.addWidget(label1)
    dialog_layout.addWidget(text1)
    dialog_layout.addWidget(label2)
    dialog_layout.addWidget(text2)
    dialog_layout.addWidget(label3)
    dialog_layout.addWidget(text3)
    dialog_layout.addWidget(label4)
    dialog_layout.addWidget(text4)
    dialog_layout.addWidget(label5)
    dialog_layout.addWidget(text5)
    dialog_layout.addWidget(label6)
    dialog_layout.addWidget(text6)
    dialog_layout.addWidget(label8)
    dialog_layout.addWidget(text8)
    dialog_layout.addWidget(label7)
    dialog_layout.addWidget(add_button)
    dialog_layout.addWidget(back_button)

    dialog.setLayout(dialog_layout)
    dialog.exec_()
     
def add_observatory_dialog(name="TUG", long="-32.779583400", lat="36.825",gmt="+3", edit=False):
    dialog = QDialog()
    dialog.closeEvent = close_dialog
    dialog_layout = QVBoxLayout()
    dialog.setWindowTitle('Add New Observatory')

    label1 = QLabel('Name')
    text1 = QLineEdit(str(name))
    label2 = QLabel("long(deg) ('-' for east)")
    text2 = QLineEdit(str(long))
    label3 = QLabel("lat(deg) ('-' for south)")
    text3 = QLineEdit(str(lat))
    label4 = QLabel('TBZ - GMT')
    text4 = QLineEdit(str(gmt))
    label7 = QLabel()
    add_button = QPushButton("add")
    add_button.clicked.connect(lambda: new_obs_add(text1.text(), text2.text(), text3.text(), text4.text(), label7, dialog))
    edit_button = QPushButton("edit")
    edit_button.clicked.connect(lambda: new_obs_edit(text1.text(), text2.text(), text3.text(), text4.text(), label7, dialog))
    back_button = QPushButton("back")
    back_button.clicked.connect(lambda: back(dialog))

    dialog_layout.addWidget(label1)
    dialog_layout.addWidget(text1)
    dialog_layout.addWidget(label2)
    dialog_layout.addWidget(text2)
    dialog_layout.addWidget(label3)
    dialog_layout.addWidget(text3)
    dialog_layout.addWidget(label4)
    dialog_layout.addWidget(text4)
    dialog_layout.addWidget(label7)
    if not edit:
        dialog_layout.addWidget(add_button)
    else:
        dialog_layout.addWidget(edit_button)
    dialog_layout.addWidget(back_button)

    dialog.setLayout(dialog_layout)
    dialog.exec_()
    
def close_dialog(close_event):
    refresh()

def new_star_add(name, ra, dec, E, P ,mag, mintime_start, label7, dialog):
    if star_file.endswith(".csv"):
        try:
            label7.setStyleSheet("color: green;")
            stars.append(Star(name, ra, dec, E, P, mag))
            label7.setText("Star Added")
            star_data = open(star_file, "a")
            star_data.write(f"{name},{ra},{dec},{E},{P},{mag},{mintime_start}\n")
            star_data.close()
            refresh()
        except Exception as e:
            label7.setStyleSheet("color: red;")
            label7.setText(str(e))

def new_obs_add(name, long, lat, tbz_gmt, label7, dialog):
    global observatories
    if observatory_file.endswith(".csv"):
        try:
            label7.setStyleSheet("color: green;")
            observatories.append(Location(name, long , lat, tbz_gmt))
            label7.setText("obs Added")
            obs_data = open(observatory_file, "a")
            obs_data.write(f"{name},{long},{lat},{tbz_gmt}\n")
            obs_data.close()
            dlg.observatory_select.addItems([f"{len(observatories)-1}: {observatories[len(observatories)-1].name}"])
            refresh()
            dialog.close()
        except Exception as e:
            label7.setStyleSheet("color: red;")
            label7.setText(str(e))
            
def new_obs_edit(name, long, lat, tbz_gmt, label7, dialog):
    global observatories
    if observatory_file.endswith(".csv"):
        try:
            label7.setStyleSheet("color: green;")
            new_obs = Location(name, long , lat, tbz_gmt)
            selected_row = obs_window.obs_table.currentRow()
            observatories[selected_row] = new_obs
            label7.setText("obs Added")
            
            with open(observatory_file, 'r') as file:
                lines = file.readlines()
            if 1 <= selected_row <= len(lines):
                lines[selected_row + 1] = f"{name},{long},{lat},{tbz_gmt}\n"
            with open(observatory_file, 'w') as file:
                file.writelines(lines)
            dlg.observatory_select.addItems([f"{observatories[len(observatories)-1].name}"])
            refresh()
            dialog.close()
            

        except Exception as e:
            label7.setStyleSheet("color: red;")
            label7.setText(str(e))
        
def back(dialog):
    dialog.close()
        
def refresh():
    global df
    global M
    global start_date
    global observer
    
    sorted_column_index = dlg.star_table.horizontalHeader().sortIndicatorSection()
    sort_order = dlg.star_table.horizontalHeader().sortIndicatorOrder()
    dlg.star_table.horizontalHeader().setSortIndicator(-1, Qt.AscendingOrder)
    
    #get the selected items row
    selected_items = dlg.star_table.selectedItems()
    if selected_items:
        selected_row = dlg.star_table.selectedItems()[0].row()
        selected_star = dlg.star_table.item(selected_row,0).text()
    else:
        selected_star = None
    
    dataframe_refresh()
    nrows, ncols = df.shape
    if not df.empty:
        df = df.sort_values(by = ["obs", "_mintime_start_h"], ascending=[False, True])
        dlg.star_table.setColumnCount(len(df.columns)-2)
        dlg.star_table.setRowCount(len(df))
        dlg.star_table.setHorizontalHeaderLabels(df.columns)
        header = dlg.star_table.horizontalHeader()  
        header.setSectionResizeMode(0, QHeaderView.Stretch)
    
    #column size
    max_len = []
    for row in df:
        l = 0
        for i in df[row].values:
            if len(str(i)) > l:
                l = len(str(i))
        if (len(row) > l):
            l = len(row)
        max_len.append(l)
    for i in range(len(df.columns)):
        dlg.star_table.setColumnWidth(i, max_len[i]*9)

    for row, row_data in enumerate(df.values):
        for col, cell_data in enumerate(row_data):
            item = QTableWidgetItem(str(df.iloc[row, col]))
            if isinstance(cell_data, (int, float)):
                item.setData(Qt.EditRole, float(cell_data))
            dlg.star_table.setItem(row, col, item)
            if df.iloc[row]["obs"]:
                item.setBackground(QColor("light green"))
            else:
                item.setBackground(QColor(255, 192, 192))
            
    
    if selected_star is not None:
        for row in range(dlg.star_table.rowCount()):
            item_text = dlg.star_table.item(row, 0).text()
            if item_text == selected_star:
                on_cell_clicked(row, 0)
                break
            
    if sort_order == 0:
        dlg.star_table.horizontalHeader().setSortIndicator(sorted_column_index, Qt.AscendingOrder)
    else:
        dlg.star_table.horizontalHeader().setSortIndicator(sorted_column_index, Qt.DescendingOrder)
        
def date_changed(value, datetype):
    global start_date
    global day_changed
    if day_changed == True:
        if datetype==0:
            if value - start_date.day == 1:
                start_date += timedelta(days=1)
                correct_date()
            elif value - start_date.day == -1:
                start_date -= timedelta(days=1)
                correct_date()
            else:
                try:
                    start_date = start_date.replace(day=value)
                except:
                    correct_date()
        if datetype==1:
            if value - start_date.month == 1:
                start_date += timedelta(days=30)
                correct_date()
            elif value - start_date.month == -1:
                start_date -= timedelta(days=30)
                correct_date()
            else:
                try:
                    start_date = start_date.replace(month=value)
                except:
                    correct_date()
        if datetype==2:
            try:
                start_date = start_date.replace(year=value)
            except:
                start_date.replace(year=value, day=start_date.day-1)
        if datetype==3:
            if value - start_date.hour == 1:
                start_date += timedelta(hours=1)
                correct_date()
            elif value - start_date.hour == -1:
                start_date -= timedelta(hours=1)
                correct_date()
            else:
                try:
                    start_date = start_date.replace(hour=value)
                except:
                    correct_date()
        if datetype==4:
            if value - start_date.minute == 1:
                start_date += timedelta(minutes=1)
                correct_date()
            elif value - start_date.minute == -1:
                start_date -= timedelta(minutes=1)
                correct_date()
            else:
                try:
                    start_date = start_date.replace(minute=value)
                except:
                    correct_date()
        refresh()
    
def correct_date():
    global day_changed
    day_changed = False
    dlg.sdd.setValue(start_date.day)
    dlg.sdm.setValue(start_date.month)
    dlg.sdy.setValue(start_date.year)
    dlg.sdhour.setValue(start_date.hour)
    dlg.sdmin.setValue(start_date.minute)
    day_changed = True
    
def moon_h_calculator(times):
    hs_moon = []
    for t in times:
        M = ephem.Moon(observer)
        M.compute(observer)
        stars[-1].ra = stars[-1]._hours_to_h(str(M.ra))
        stars[-1].dec = stars[-1]._deg_to_float(str(M.dec))
        hs_moon.append(h_calculator(t, stars[-1]))
    return hs_moon
        
def time_h_plot(star_name):
    global s
    global min_alt
    global twilight_start
    global twilight_end
    
    times = []
    
    delta = timedelta(minutes=1)
    t = twilight_start - timedelta(hours=2)
    while t <= twilight_end + timedelta(hours=2):
        times.append(t)
        t += delta
        
    # plotting star find in list
    for i in stars:
        if i.name == star_name:
            plotting_star = i
            break
    else:
        ValueError("Couldn't find star")
        
    fig, ax = plt.subplots()
    hs = [h_calculator(time, plotting_star) for time in times]
    hs_moon = moon_h_calculator(times)
    
    #plot hours
    ax.plot(times, hs, "r-", label="Star")
    ax.plot(times, hs_moon, "k--", linewidth=1, label="Moon")
    ax.axhline(y=min_alt, color='b', linestyle='-', label="Min Height")

    # min times axvlines
    for i in plotting_star.mintimes:
        mintime1 = np.datetime64(i)
        ax.axvline(x=mintime1, color="g", alpha=.4)
    for i in plotting_star.mintimes2:
        mintime2 = np.datetime64(i)
        ax.axvline(x=mintime2, color="y", alpha=.4)
    all_min_times = np.append(plotting_star.mintimes, plotting_star.mintimes2)
    ax2 = ax.secondary_xaxis("top")
    ax2.tick_params(axis="x", color="red")
    ax2.set_xticks([*all_min_times], minor=False)
    min_labels = [f"{t.hour:02d}:{t.minute:02d}" for t in all_min_times]
    ax2.set_xticklabels(min_labels)
    ax2.tick_params(labelsize=8)
        
    ax2.set_xlabel("min times")

    # label namimg
    ax.set_ylabel("h(deg)")
    ax.set_xlabel("hour")
    
    ax.axvline(twilight_start, color="m", linestyle="--")
    ax.axvline(twilight_end, color="m", linestyle="--")
    
    #format x axis to hour min
    bx = plt.gca()
    bx.xaxis.set_major_locator(mdates.AutoDateLocator())
    bx.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    bx.tick_params(axis='x', labelsize=8)
    
    ax.set_ylim(0, 90)
    ax.set_xlim(twilight_start - timedelta(minutes=30), twilight_end + timedelta(minutes=30))
    
    for i in plotting_star.mintimes:
        if ~np.isnan(plotting_star.minstart):
            print(plotting_star.minstart)
            som = i - timedelta(minutes=plotting_star.minstart)
            eom = i + timedelta(minutes=plotting_star.minstart)
        else:
            som = i - timedelta(minutes=mst)
            eom = i + timedelta(minutes=mst)
        x = som + np.arange(10) * (eom - som) / (10 - 1)
        y_min = -200
        y_max = 200
        plt.fill_between(x, y_min, y_max, color='gray', alpha=0.3)
        
    current_h = h_calculator(start_date, plotting_star)
    ax.plot(start_date, current_h, 'go', label="Current Time") 
    dlg.star_name.setText(star_name)
    
    #legend
    if legend:
        mintime1_legend = mlines.Line2D([], [], color='g', alpha=.4, label="Mintimes 1", linestyle='-')
        mintime2_legend = mlines.Line2D([], [], color='y', alpha=.4, label="Mintimes 2", linestyle='-')
        night_legend = mlines.Line2D([], [], color='m', alpha=.4, label="Night Start/End", linestyle='--')
        plt.legend()
        handles, labels = ax.get_legend_handles_labels()
        handles.append(mintime1_legend)
        handles.append(mintime2_legend)
        handles.append(night_legend)
        plt.legend(handles=handles)

    if dlg.graph_layout.count() == 0:
        canvas = FigureCanvas(fig)
        dlg.graph_layout.addWidget(canvas)
        toolbar = NavigationToolbar(canvas, dlg.graph_grid)
        dlg.graph_layout.addWidget(toolbar)
    else:
        for i in reversed(range(dlg.graph_layout.count())):
            widget = dlg.graph_layout.itemAt(i).widget()
            if widget is not None:
                widget.close()
            dlg.graph_layout.takeAt(i)

        canvas = FigureCanvas(fig)
        dlg.graph_layout.addWidget(canvas)
        toolbar = NavigationToolbar(canvas, dlg.graph_grid)
        dlg.graph_layout.addWidget(toolbar)

    plt.close()

def on_cell_clicked(row, column):
    global lightcurve
    global monthing_star
    dlg.star_table.selectRow(row)
    selected_row = dlg.star_table.currentRow()
    item = dlg.star_table.item(selected_row, 0)
    monthing_star = item.text()
    time_h_plot(monthing_star)
    if lightcurve:
        try:
            worker = Worker(lambda: tess_analysis(monthing_star))

            # Connect a slot to the worker's finished signal
            worker.signals.finished.connect(plot_lightcurve)

            # Start the worker in the global QThreadPool
            pool = QThreadPool.globalInstance()
            pool.start(worker)
        except:
            pass
    

def plot_lightcurve(result):
    global monthing_star
    if not result:
        dlg.async_label.setText(f"{monthing_star} has no TESS lightcurve")
        return
    tess_times = result[0]
    tess_sap_flux = result[1]
    tess_sap_flux_err = result[2]
    tess_t = result[3]
    plotting_star_tess = result[4]
    if monthing_star == plotting_star_tess.name and int(dlg.graph_layout.count()) < 4:
        fig, ax = plt.subplots()
        ax.plot(tess_times, tess_sap_flux, "r.")
        ax.set_xlim(min(tess_t), max(tess_t))
        canvas = FigureCanvas(fig)
        dlg.graph_layout.addWidget(canvas)
        toolbar = NavigationToolbar(canvas)
        dlg.graph_layout.addWidget(toolbar)
        dlg.async_label.setText(f"{monthing_star} TESS lightcurve downloded")
        plt.close()
    
def observatory_selected(index):
    global obs_loc
    selected_text = dlg.observatory_select.currentIndex()
    obs_loc = observatories[selected_text]
    if(str(obs_loc.tbz_gmt)[0] != "-"):
        dlg.gmt_text.setText(f"GMT: +{obs_loc.tbz_gmt}")
    refresh()
    
def msp_changed():
    global mst
    text = dlg.min_start_period_input.text()
    try:
        msp = float(text)
        if msp >= 1000:
            raise ValueError("msp can't be more than 1")
        elif msp < 0:
            raise ValueError("msp can't be less than 0")
        else:
            mst = msp
        settings_changed()
        refresh()
    except:
        dlg.min_start_period_input .setText("60")

def min_alt_changed():
    global min_alt
    try:
        ma = float(dlg.min_altitude.text())
        if ma < 0:
            raise ValueError("min altitude can't be less than 0")
        elif ma > 90:
            raise ValueError("min altitude can't be more than 90")
        else:
            min_alt = ma
            refresh()
        settings_changed()
    except:
        dlg.min_altitude.setText("30")
        
def deselect_button_clicked():
    dlg.star_table.clearSelection()
    while dlg.graph_layout.count():
        item = dlg.graph_layout.itemAt(0)
        widget = item.widget()
        dlg.graph_layout.removeItem(item)
        widget.setParent(None)

def go_time(given_time):
    global start_date
    if given_time == "now":
        start_date = datetime.now().replace(second=0)
        correct_date()
        refresh()
        return
    elif given_time == "tonight":
        if datetime.now().hour < 8:
            desired_date = date.today() - timedelta(days=1)
        else:
            desired_date = date.today()

    elif given_time == "selected night":
        if start_date.hour < 8:
            desired_date = start_date - timedelta(days=1)
        else:
            desired_date = start_date
            
    desired_time = time(20, 0, 0)
    start_date = datetime.combine(desired_date, desired_time)
    correct_date()
    refresh()
    
def legend_box_clicked(state):
    global legend
    if state == 2: 
        legend = True
        refresh()
    else:
        legend = False
        refresh()
        
def lightcurve_box_clicked(state):
    global lightcurve
    if state == 2: 
        lightcurve = True
        refresh()
    else:
        lightcurve = False
        refresh()
        
def sector_query(star_name):
    try:
        query = Tesscut.get_sectors(objectname=star_name)
        sectors = [sector for sector in query["sector"]]
        return sectors  
    except:
        pass
        
def tess_analysis(object_name):
    dlg.async_label.setText(f"Downloading {object_name} lightcurve")
    folder_path = "lightcurves"
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    else:
        pass
    file_name = f"{object_name}.txt"
    files = [name for name in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, name))]
    
    for i in stars:
        if i.name == object_name:
            plotting_star_tess = i
            break
    else:
        dlg.async_label.setText("Couldn't find lightcurve")

    if file_name not in files:
        sector_numbers = sector_query(object_name)
        for sn in sector_numbers[::-1]:
            print(sn)
            try:
                sector = sn
                query = Observations.query_criteria(obs_collection="TESS",
                                                    objectname=object_name,
                                                    dataproduct_type="timeseries",
                                                    sequence_number=sector,
                                                    t_exptime = 120,
                                                    radius=0)

                products = Observations.get_product_list(query)
                filtered = Observations.filter_products(products,
                                                        productSubGroupDescription="LC")

                downloaded_data = Observations.download_products(filtered)
                path = downloaded_data["Local Path"][0]
                hdul = fits.open(path)
                file = open(folder_path + "/" + file_name, "w")
                data = ""
                """
                flux = flux[mask]
                error = error[mask]
                """
                for j in range(len(hdul[1].data["TIME"])):
                    data += str(hdul[1].data["TIME"][j]) + " " + str(hdul[1].data["PDCSAP_FLUX"][j]) + " " + str(hdul[1].data["PDCSAP_FLUX_ERR"][j]) + "\n"
                file.write(data)
                file.close()
                hdul.close()
                os.remove(path)
                break
            except Exception as e:
                pass
        else:
            file = open(folder_path + "/" + file_name, "w")
            file.write("No Tess Found")
        
    try:
        file = np.loadtxt(folder_path + "/" + file_name)
        tess_times = file[:, 0]
        tess_sap_flux = file[:, 1]
        tess_sap_flux_err = file[:, 2]
    except:
        return False
    
    tess_t = tess_times
    
    mid = int(len(tess_t)/2)
    start_cand = tess_t[mid]
    i = mid
    a = start_cand + plotting_star_tess.P * 2
    while tess_t[i] < start_cand + plotting_star_tess.P * 2 or np.isnan(tess_t[i]):
        if i == len(tess_t)-1:
            i = tess_t[mid]
            break
        if np.isnan(tess_t[i]) or np.isnan(tess_sap_flux[i]):
            i = i+1
            while np.isnan(tess_t[i]):
                i = i+1
            start_cand = tess_t[i]
        else:
            i = i+1
    else:
        start_time = start_cand
        
    mask = np.where((tess_t > start_time) & (tess_t < start_time + plotting_star_tess.P * 2))
    tess_t = tess_t[mask]
    
    return tess_times, tess_sap_flux, tess_sap_flux_err, tess_t, plotting_star_tess

def on_cell_double_clicked(item):
    global star_for_month
    found = 1
    month_min_altitude.setText(str(min_alt)) 
    
    row = dlg.star_table.selectedItems()[0].row()
    star_name = dlg.star_table.item(row, 0).text()
    for i in stars:
        if i.name == star_name:
            star_for_month = i
            break
    else:
        found = 0
    if found == 1:
        month_refresh()
        month.show()
        month.activateWindow()
    
def month_refresh():
    global star_for_month
    monthing_star = star_for_month
    month.month_star_name.setText(monthing_star.name)
    month.month_ra.setText(str(monthing_star.ra_str))
    month.month_dec.setText(str(monthing_star.dec_str))
    month_mag.setText(str(monthing_star.mag))
    month.month_phase.setText(str(monthing_star.P))
    
        
    fred = ephem.Observer()
    if start_date.hour < 12:
        fred.date = start_date.replace(hour=12) - timedelta(days=1)
    else:
        fred.date = start_date.replace(hour=12)

    fred.lon = str(obs_loc.long_h*15*-1) 
    fred.lat = str(obs_loc.lat_deg)      
    fred.elev = 20
    fred.pressure= 0
    fred.horizon = '-18'
    
    #start_end_date
    datess = []
    phase_df = pd.DataFrame()
    alt_df = pd.DataFrame()
    s_date = datetime(month_start_date.date().year(), month_start_date.date().month(), month_start_date.date().day(), month_start_date.time().hour(), minute=month_start_date.time().minute())
    e_date = datetime(month_end_date.date().year(), month_end_date.date().month(), month_end_date.date().day(), month_end_date.time().hour(), minute=month_end_date.time().minute())
    
    date_delta = timedelta(days = 1)
    d = s_date
    while d <= e_date:
        datess.append(d)
        d += date_delta
    earliest_tw_start = 23
    latest_tw_end = 1
    for i in datess:
        fred.date = i.replace(hour=12)
        tw_end=fred.next_rising(ephem.Sun(), use_center=True).datetime() + timedelta(hours=obs_loc.tbz_gmt)
        tw_start=fred.next_setting(ephem.Sun(), use_center=True).datetime() + timedelta(hours=obs_loc.tbz_gmt)
        if tw_start.hour < earliest_tw_start:
            earliest_tw_start = tw_start.hour
        if tw_end.hour > latest_tw_end:
            latest_tw_end = tw_end.hour
    for i in datess:
        t = i.replace(hour=earliest_tw_start)
        while t <= (i + timedelta(days=1)).replace(hour=latest_tw_end):
            t += timedelta(hours=1)
            phase_df.loc[f"{i.day:02d}.{i.month:02d}.{i.year:02d}", f"{t.hour:02d}:00"] = phase_calculator(t, monthing_star)
            alt_df.loc[f"{i.day:02d}.{i.month:02d}.{i.year:02d}", f"{t.hour:02d}:00"] = h_calculator(t, monthing_star)
    
    month.month_table_widget.setRowCount(phase_df.shape[0])
    month.month_table_widget.setColumnCount(phase_df.shape[1])
    month.month_table_widget.setHorizontalHeaderLabels(phase_df.columns)
    month.month_table_widget.setVerticalHeaderLabels(phase_df.index)

    for row in range(phase_df.shape[0]):
        mintimes = []
        d, m, y = phase_df.index[row].split(".")
        time = s_date.replace(day=int(d), month=int(m), year=int(y))
        gst, jd = gst_jd_calculator(time)
        ph = ((jd-monthing_star.E)/monthing_star.P)%1
        first = time.replace(hour=earliest_tw_start)
        last = time.replace(hour=latest_tw_end) + timedelta(days=1)
        while time > first:
            time -= timedelta(hours=(ph)*monthing_star.P*24)
            ph = 1
            if time > first and time > first:
                mintimes.append(time)
        ph = ((jd-monthing_star.E)/monthing_star.P)%1
        time = s_date.replace(day=int(d), month=int(m), year=int(y))
        while time < last:
            time += timedelta(hours=(1-ph)*monthing_star.P*24)
            ph = 0
            if time < last and time > first:
                mintimes.append(time)
                
        mintimes2 = []
        d, m, y = phase_df.index[row].split(".")
        time = s_date.replace(day=int(d), month=int(m), year=int(y))
        gst, jd = gst_jd_calculator(time)
        ph = ((jd-monthing_star.E)/monthing_star.P)%1 + .5
        while time > first:
            time -= timedelta(hours=(ph)*monthing_star.P*24)
            ph = 1
            if time > first and time > first:
                mintimes2.append(time)
        ph = ((jd-monthing_star.E)/monthing_star.P)%1 + .5
        time = s_date.replace(day=int(d), month=int(m), year=int(y))
        while time < last:
            time += timedelta(hours=(1-ph)*monthing_star.P*24)
            ph = 0
            if time < last and time > first:
                mintimes2.append(time)
        
        mintimes_hour = []
        mintimes_hour2 = []
        for i in mintimes:
            mintimes_hour.append(find_closest_hour(i))
        for i in mintimes2:
            mintimes_hour2.append(find_closest_hour(i))
        
        for col in range(phase_df.shape[1]):
            phase = phase_df.iloc[row][col]
            alt = alt_df.iloc[row][col]
            item = QTableWidgetItem(f"{phase_df.iloc[row, col]:.4f}")
            month.month_table_widget.setItem(row, col, item)
            if float(alt) > month_min_alt:
                item.setBackground(QColor(150, 250, 150))
            else:
                item.setBackground(QColor(250, 150, 150))
            
            if int(phase_df.columns[col][:2]) in mintimes_hour:
                closest_item = month.month_table_widget.item(row, col)
                closest_item.setForeground(QColor(60, 60, 250))
            if int(phase_df.columns[col][:2]) in mintimes_hour2:
                closest_item = month.month_table_widget.item(row, col)
                closest_item.setForeground(QColor(150, 0, 250))
                
    
    for row_index in range(month.month_table_widget.rowCount()):
        set_row_alignment(month.month_table_widget, row_index, Qt.AlignCenter)

def find_closest_hour(dt):
    hour = dt.hour
    minutes = dt.minute

    if minutes >= 30:
        closest_hour = (hour + 1) % 24
    else:
        closest_hour = hour

    return closest_hour
    
def set_row_alignment(table_widget, row_index, alignment):
    for col in range(table_widget.columnCount()):
        item = table_widget.item(row_index, col)
        if item:
            item.setTextAlignment(alignment)
            
def phase_calculator(gdate, star):
    gst, jd = gst_jd_calculator(gdate)
    ph = ((jd-star.E)/star.P)%1
    return ph     

def month_min_altitude_changed(value):
    global month_min_alt
    try:
        value = int(value)
        month_min_alt = value
        month_refresh()
        settings_changed()
    except:
        month_min_altitude.setText(str(min_alt))     
        
def on_month_double_clicked(item):
    global star_for_month
    global start_date
    month_min_altitude.setText(str(min_alt)) 
    row = item.row()
    col = item.column()
    row_label = month.month_table_widget.verticalHeaderItem(row).text()
    column_label = month.month_table_widget.horizontalHeaderItem(col).text()
    d, m, y = row_label.split(".")
    h, min = column_label.split(":")
    date = datetime(int(y), int(m), int(d), int(h), int(min))
    if int(h)<12:
        date = date + timedelta(days=1)
    dt = date
    start_date = dt
    correct_date()
    for row in range(dlg.star_table.rowCount()):
        item = dlg.star_table.item(row, 0) 
        if item.text() == star_for_month.name:
            on_cell_clicked(row,0)
    refresh()
    dlg.show()
    dlg.activateWindow()
    
def edit_star():
    global star_for_month
    dialog = QDialog()
    dialog_layout = QVBoxLayout()
    dialog.setWindowTitle('Add New Star')
    dialog.closeEvent = close_dialog
    label1 = QLabel('Name')
    text1 = QLineEdit(star_for_month.name)
    text1.setDisabled(True)
    label2 = QLabel('RA')
    text2 = QLineEdit(str(month.month_ra.text()))
    label3 = QLabel('Dec')
    text3 = QLineEdit(str(month.month_dec.text()))
    label4 = QLabel('E')
    text4 = QLineEdit(str(star_for_month.E))
    label5 = QLabel('P')
    text5 = QLineEdit(str(star_for_month.P))
    label6 = QLabel('M')
    text6 = QLineEdit(str(star_for_month.mag))
    label8 = QLabel('Mintime - Mintime Start (minute)(optional)')
    text8 = QLineEdit(str(star_for_month.minstart))
    label7 = QLabel()
    add_button = QPushButton("edit")
    add_button.clicked.connect(lambda: star_edited(text1.text(), text2.text(), text3.text(), text4.text(), text5.text(), text6.text(), text8.text(), label7))
    back_button = QPushButton("back")
    back_button.clicked.connect(lambda: back(dialog))
    
    dialog_layout.addWidget(label1)
    dialog_layout.addWidget(text1)
    dialog_layout.addWidget(label2)
    dialog_layout.addWidget(text2)
    dialog_layout.addWidget(label3)
    dialog_layout.addWidget(text3)
    dialog_layout.addWidget(label4)
    dialog_layout.addWidget(text4)
    dialog_layout.addWidget(label5)
    dialog_layout.addWidget(text5)
    dialog_layout.addWidget(label6)
    dialog_layout.addWidget(text6)
    dialog_layout.addWidget(label8)
    dialog_layout.addWidget(text8)
    dialog_layout.addWidget(label7)
    dialog_layout.addWidget(add_button)
    dialog_layout.addWidget(back_button)

    dialog.setLayout(dialog_layout)
    dialog.exec_()
    
def star_edited(name, ra, dec, E, P ,mag, minstart, label7):
    if star_file.endswith(".csv"):
        try:
            label7.setStyleSheet("color: green;")
            old_name = star_for_month.name
            star_for_month.name = name
            star_for_month.ra = star_for_month._hours_to_h(ra)
            star_for_month.dec = star_for_month._deg_to_float(dec)
            star_for_month.E = float(E)
            star_for_month.P = float(P)
            star_for_month.mag = mag
            star_for_month.minstart = minstart
            label7.setText("Star Edited")
            
            with open(star_file, 'r') as file:
                lines = file.readlines()
            for i, line in enumerate(lines):
                if line[:len(old_name)] == old_name:
                    target_line_number = i
                    break
            file.close()
            
            if target_line_number >= 1 and target_line_number <= len(lines):
                lines[target_line_number] = f"{name},{ra},{dec},{E},{P},{mag},{minstart}\n"
                with open(star_file, 'w') as file:
                    file.writelines(lines)
            
            stars = []
            dlg.star_table.clearContents()
            load()
            add_moon()
            refresh()
        except Exception as e:
            label7.setStyleSheet("color: red;")
            label7.setText(str(e))
        
def delete_star():
    global stars
    global star_for_month
    global df
    with open(star_file, 'r') as file:
        lines = file.readlines()
    for i, line in enumerate(lines):
        if line[:len(star_for_month.name)] == star_for_month.name:
            target_line_number = i
            break
    file.close()
    
    if target_line_number >= 1 and target_line_number <= len(lines):
        del lines[target_line_number]

        with open(star_file, 'w') as file:
            file.writelines(lines)
    
    month.close()  
    stars = []
    dlg.star_table.clearContents()
    df = pandas.DataFrame()
    load()
    add_moon()
    refresh()
    
def add_moon():
    global observer
    observer = ephem.Observer()
    observer.lat = obs_loc.lat_deg  # Latitude of the observer (e.g., London)
    observer.lon = obs_loc.long_h
    observer.date = f"{start_date.year}/{start_date.month}/{start_date.day} {start_date.hour}:{start_date.minute}"
    M = ephem.Moon()
    M.compute(observer)
    stars.append(Star("moon", str(M.ra), str(M.dec), E=60099.65417, P=29.530588, mag=0, pmra=0, pmdec=0))
  
def go_to_simbad():
    try:
        chrome_options = Options()
        chrome_options.add_argument('--headless') 

        driver = webdriver.Chrome(options=chrome_options) 

        driver.get('https://simbad.cds.unistra.fr/simbad/sim-fbasic')
        search_input = driver.find_element(By.NAME, 'Ident')
        search_input.send_keys(monthing_star)
        search_form = driver.find_element(By.NAME, 'submit')
        search_form.click()

        current_url = driver.current_url
        driver.quit()
        
        webbrowser.open_new_tab(current_url)
    except:
        pass
    
def load_stars():
    global star_file
    star_file_old = star_file
    options = QFileDialog.Options()
    options |= QFileDialog.DontUseNativeDialog 
    
    file_dialog = QFileDialog()
    file_dialog.setOptions(options)
    
    file_dialog.setFileMode(QFileDialog.ExistingFile)
    
    if file_dialog.exec_():
        selected_files = file_dialog.selectedFiles()
        star_file = selected_files[0]
        
        try:
            load()
            deselect_button_clicked()
        except Exception as e:
            print(e)
            message_box = QMessageBox()
            message_box.setIcon(QMessageBox.Critical)
            message_box.setWindowTitle("Error")
            message_box.setText(str(e))
            message_box.exec_()
            star_file = star_file_old
            load()
        settings_changed()
        add_moon()
        refresh()
    
def load_observatories():
    global observatory_file
    observatory_file_old = observatory_file
    options = QFileDialog.Options()
    options |= QFileDialog.DontUseNativeDialog 
    
    file_dialog = QFileDialog()
    file_dialog.setOptions(options)
    
    file_dialog.setFileMode(QFileDialog.ExistingFile)
    
    if file_dialog.exec_():
        selected_files = file_dialog.selectedFiles()
        observatory_file = selected_files[0]
        try:
            load()
        except:
            message_box = QMessageBox()
            message_box.setIcon(QMessageBox.Critical)
            message_box.setWindowTitle("Error")
            message_box.setText("Observatory file couldn't be loaded")
            message_box.exec_()
            observatory_file = observatory_file_old
            load()
        settings_changed()
        add_moon()
        refresh()
        
def settings_changed():
    global columns_active
    if os.path.exists("settings.dat"):
        os.remove("settings.dat")
    settings = {
    'star_file': star_file,
    'observatory_file': observatory_file,
    'mst': mst,
    'min_alt': min_alt,
    'columns_active':columns_active
    }
    with open('settings.dat', 'wb') as file:
        pickle.dump(settings, file)
        
def sorting_changed():

    sorting = dlg.sorting_combobox.currentText()
    if sorting == "Sorting: Default":
        dlg.star_table.setSortingEnabled(False)
    else:
        dlg.star_table.setSortingEnabled(True)
    refresh()
    
def on_item_selection_changed():
    selected_items = dlg.star_table.selectedItems()
    if selected_items:
        dlg.clear_button.setDisabled(False)
        dlg.star_page_button.setDisabled(False)
    else:
        dlg.clear_button.setDisabled(True)
        dlg.star_page_button.setDisabled(True)

def columns_clicked():
    global columns_active
    global columns_active_old
    global col_dialog
    columns_active_old = columns_active
    col_dialog = QDialog()
    col_dialog.setWindowTitle('Checkbox Dialog')

    layout = QVBoxLayout()

    checkboxes = []
    
    current_alt_checkbox = QCheckBox('Current Altitude', col_dialog)
    layout.addWidget(current_alt_checkbox)
    checkboxes.append(current_alt_checkbox)
    current_phase_checkbox = QCheckBox('Current Phase', col_dialog)
    layout.addWidget(current_phase_checkbox)
    checkboxes.append(current_phase_checkbox)
    period_checkbox = QCheckBox('Period', col_dialog)
    layout.addWidget(period_checkbox)
    checkboxes.append(period_checkbox)
    min_time_start_checkbox = QCheckBox('Min Time Start', col_dialog)
    layout.addWidget(min_time_start_checkbox)
    checkboxes.append(min_time_start_checkbox)
    min_time_end_checkbox = QCheckBox('Min Time End', col_dialog)
    layout.addWidget(min_time_end_checkbox)
    checkboxes.append(min_time_end_checkbox)
    min_time_checkbox = QCheckBox('Min Time', col_dialog)
    layout.addWidget(min_time_checkbox)
    checkboxes.append(min_time_checkbox)
    mag_checkbox = QCheckBox('Magnitude', col_dialog)
    layout.addWidget(mag_checkbox)
    checkboxes.append(mag_checkbox)
    mintime_h_checkbox = QCheckBox('Altitude At Mintime', col_dialog)
    layout.addWidget(mintime_h_checkbox)
    checkboxes.append(mintime_h_checkbox)
    moon_distance_checkbox = QCheckBox('Moon Distance', col_dialog)
    layout.addWidget(moon_distance_checkbox)
    checkboxes.append(moon_distance_checkbox)
    
    for i, checkbox in enumerate(checkboxes):
        checkbox.setChecked(columns_active[i])
        checkbox.stateChanged.connect(lambda state, checkbox=checkbox, i=i: on_checkbox_state_changed(checkbox, columns_active, i))

    save_button = QPushButton('Save', col_dialog)
    cancel_button = QPushButton('Cancel', col_dialog)

    save_button.clicked.connect(columns_save_clicked)
    cancel_button.clicked.connect(columns_cancel_clicked)

    layout.addWidget(save_button)
    layout.addWidget(cancel_button)

    col_dialog.setLayout(layout)

    result = col_dialog.exec_()
        
def on_checkbox_state_changed(checkbox, columns_active, checkbox_index):
    columns_active[checkbox_index] = checkbox.isChecked()

def columns_save_clicked():
    col_dialog.accept()
    settings_changed()
    refresh()

def columns_cancel_clicked():
    global columns_active_old
    global columns_active
    columns_active = columns_active_old
    col_dialog.reject()
    settings_changed()
    refresh()
    
def open_obs():
    obs_window.obs_table.clearContents()
    obs_window.show()
    obs_window.obs_table.setRowCount(len(observatories))
    obs_window.obs_table.setColumnCount(2)
    obs_window.obs_table.setHorizontalHeaderLabels(["Observatory", "GMT"])
    for row, location in enumerate(observatories):
        name_item = QTableWidgetItem(location.name)
        if str(location.tbz_gmt)[0] != "-":
            tbz_item = QTableWidgetItem("+" + str(location.tbz_gmt))
        else:
            tbz_item = QTableWidgetItem(str(location.tbz_gmt))
        name_item.setTextAlignment(Qt.AlignCenter)
        tbz_item.setTextAlignment(Qt.AlignCenter)
        obs_window.obs_table.setItem(row, 0, name_item)
        obs_window.obs_table.setItem(row, 1, tbz_item)
    obs_window.activateWindow()
    
def on_cell_clicked_obs(row, column):
    obs_window.obs_table.selectRow(row)
    selected_row = obs_window.obs_table.currentRow()
    obs_window.obs_lon.setText(str(observatories[selected_row].long))
    obs_window.obs_lat.setText(str(observatories[selected_row].lat_deg))
    obs_window.obs_gmt.setText("+" + str(observatories[selected_row].tbz_gmt))
    
def select_obs():
    global obs_loc
    selected_row = obs_window.obs_table.currentRow()
    dlg.observatory_select.setCurrentIndex(selected_row)
    
def edit_obs():
    global obs_loc
    selected_row = obs_window.obs_table.currentRow()
    sel_obs = observatories[selected_row]
    add_observatory_dialog(name=sel_obs.name, long=sel_obs.long_h, lat=sel_obs.lat_deg, gmt=sel_obs.tbz_gmt, edit=True)
    dlg.observatory_select.setItemText(selected_row, observatories[selected_row].name)
    open_obs()
    
def del_obs():
    global observatories
    selected_row = obs_window.obs_table.currentRow()
    observatories = observatories[:selected_row] + observatories[selected_row + 1:]
    obs_window.obs_table.removeRow(selected_row)

    with open(observatory_file, 'r') as file:
        lines = file.readlines()
    if 1 <= selected_row <= len(lines):
        lines.pop(selected_row + 1)
    with open(observatory_file, 'w') as file:
        file.writelines(lines)
        
    dlg.observatory_select.removeItem(selected_row)
        
    refresh()
    
def add_obs():
    global obs_loc
    add_observatory_dialog()
    open_obs()
    
def about():
    about_dialog = AboutDialog()
    about_dialog.exec_()

    
    
if __name__ == "__main__":
    app = QtWidgets.QApplication([])
    dlg = uic.loadUi("prog.ui")
    dlg.star_table.clearContents()
    
    #gloabal variables
    if os.path.exists("settings.dat"):
        with open('settings.dat', 'rb') as file:
            loaded_settings = pickle.load(file)
        if os.path.isfile(loaded_settings['star_file'].split('/')[-1]):
            star_file = loaded_settings['star_file'].split('/')[-1]
        else:
            star_file = loaded_settings['star_file']
        if os.path.isfile(loaded_settings['observatory_file'].split('/')[-1]):
            observatory_file = loaded_settings['observatory_file'].split('/')[-1]
        else:
            observatory_file = loaded_settings['observatory_file']
        mst = float(loaded_settings['mst'])
        min_alt = float(loaded_settings['min_alt'])
        columns_active = list(loaded_settings['columns_active'])
    else:
        settings = {
        'star_file': 'star.csv',
        'observatory_file': 'obs.csv',
        'mst': 60,
        'min_alt': 30,
        'columns_active': [False, False, True, True, True, True, True, True, True]
        }
        star_file = "star.csv"
        observatory_file = "obs.csv"
        mst = 60
        min_alt = 30
        columns_active = [False, False, True, True, True, True, True, True, True]
        with open('settings.dat', 'wb') as file:
            pickle.dump(settings, file)
            
    observatories = []
    stars = []
    obs_loc = None
    #mintime start period
    df = pandas.DataFrame()
    start_date = datetime.now()
    times = []
    s = None
    month_min_alt = min_alt
    
    try:
        load()
    except:
        print("settings couldn't loaded returning default settings")
        os.remove("settings.dat")
        star_file = "star.csv"
        observatory_file = "obs.csv"
        if not os.path.isfile(star_file):
            with open(star_file, mode='w', newline='') as csv_file:
                fieldnames = ['Star Name', 'ra(h)', 'dec(deg)', 'Epoch', 'Period', 'Magnitude', 'Mintime_Start']
                writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
                writer.writeheader()
        load()
    #moon
    if(str(obs_loc.tbz_gmt)[0] != "-"):
        dlg.gmt_text.setText(f"GMT: +{obs_loc.tbz_gmt}")
    add_moon()

    #gui variables
    
    dlg.star_table.cellClicked.connect(on_cell_clicked)
    dlg.star_table.itemDoubleClicked.connect(on_cell_double_clicked)
    dlg.star_table.itemSelectionChanged.connect(on_item_selection_changed)
    dlg.tonight_button.clicked.connect(lambda: go_time("tonight"))
    dlg.now_button.clicked.connect(lambda: go_time("now"))
    dlg.selectednight_button.clicked.connect(lambda: go_time("selected night"))
    model = QStandardItemModel()
    dlg.observatory_select.setEditable(True)
    dlg.observatory_select.lineEdit().setAlignment(Qt.AlignCenter)
    dlg.observatory_select.lineEdit().setReadOnly(True)
    for observatory in observatories:
        item = QStandardItem(observatory.name)
        item.setTextAlignment(Qt.AlignCenter)  # Center-align the text
        model.appendRow(item)
    dlg.observatory_select.setModel(model)
    dlg.min_start_period_input.setText(str(mst))
    dlg.min_start_period_input.editingFinished.connect(msp_changed)
    dlg.min_altitude.setText(str(min_alt))
    dlg.min_altitude.editingFinished.connect(min_alt_changed)
    dlg.observatory_select.setStyleSheet("QComboBox { text-align: center; }")
    dlg.observatory_select.currentIndexChanged.connect(observatory_selected)
    dlg.sdd.setValue(start_date.day)
    dlg.sdm.setValue(start_date.month)
    dlg.sdy.setValue(start_date.year)
    dlg.sdhour.setValue(start_date.hour)
    dlg.sdmin.setValue(start_date.minute)
    day_changed = True
    dlg.sdd.valueChanged.connect(lambda value: date_changed(value, 0))
    dlg.sdm.valueChanged.connect(lambda value: date_changed(value, 1))
    dlg.sdy.valueChanged.connect(lambda value: date_changed(value, 2))
    dlg.sdhour.valueChanged.connect(lambda value: date_changed(value, 3))
    dlg.sdmin.valueChanged.connect(lambda value: date_changed(value, 4))
    dlg.clear_button.clicked.connect(deselect_button_clicked)
    legend = False
    dlg.legend_box.stateChanged.connect(legend_box_clicked)
    lightcurve = False
    dlg.lightcurve_box.stateChanged.connect(lightcurve_box_clicked)
    dlg.action_load_stars.triggered.connect(load_stars)
    dlg.action_load_observatories.triggered.connect(load_observatories)
    dlg.action_about.triggered.connect(about)
    dlg.sorting_combobox.addItems(["Sorting: Default", "Sorting: Free"])
    dlg.sorting_combobox.setCurrentIndex(0)
    dlg.sorting_combobox.currentIndexChanged.connect(sorting_changed)
    dlg.sorting_combobox.setStyleSheet("QComboBox { text-align: center; }")
    dlg.star_page_button.clicked.connect(lambda: on_cell_double_clicked(dlg.star_table.selectedItems()[0]))
    dlg.columns_button.clicked.connect(columns_clicked)
    dlg.add_star_button.clicked.connect(add_star_dialog)
    
    month = uic.loadUi("month.ui")
    flags = month.windowFlags()
    flags &= ~Qt.WindowMinimizeButtonHint
    month.setWindowFlags(flags)
    
    month.month_table_widget.itemDoubleClicked.connect(on_month_double_clicked)
    month_mag = month.month_mag
    month_start_date = month.month_start_date
    month_end_date = month.month_end_date
    month_start_date.setDate(start_date)
    month_end_date.setDate(month_start_date.date().addDays(90))
    month_submit_date = month.month_submit_date
    month_submit_date.clicked.connect(month_refresh)
    month_min_altitude = month.month_min_altitude
    month_min_altitude.textChanged.connect(month_min_altitude_changed)
    month_edit_star = month.edit_star
    month_delete_star = month.delete_star  
    month_edit_star.clicked.connect(edit_star)
    month_delete_star.clicked.connect(delete_star)
    month_simbad = month.month_simbad
    month_simbad.clicked.connect(go_to_simbad)
    
    month_log = month.month_log
    month_log.clicked.connect(lambda: log.open_log(month.month_star_name.text()))
    
    obs_window = uic.loadUi("obs.ui")
    dlg.observatories_button.clicked.connect(open_obs)
    obs_window.obs_table.cellClicked.connect(on_cell_clicked_obs)
    obs_window.obs_select_button.clicked.connect(select_obs) 
    obs_window.obs_edit_button.clicked.connect(edit_obs)   
    obs_window.obs_del_button.clicked.connect(del_obs)
    obs_window.obs_add_button.clicked.connect(add_obs)
    obs_window.obs_back_button.clicked.connect(obs_window.close)
    cell_selected = False
    refresh()
    
    #magnitude finder
    # from astroquery.simbad import Simbad
    # for star in stars[:-1]:

    #     custom_simbad = Simbad()
    #     custom_simbad.add_votable_fields('flux(V)', 'flux(G)')  # Magnitude in V and g bands

    #     # Query Simbad for the star
    #     result_table = custom_simbad.query_object(star.name)

    #     # Retrieve the magnitude
    #     if 'FLUX_V' in result_table.colnames:
    #         magnitude = result_table['FLUX_V'][0]
    #         band = 'V'
    #     if not isinstance(magnitude, float):
    #         magnitude = result_table['FLUX_G'][0]
    #         band = 'g'
    #     print(magnitude)

    
    #gui start
    dlg.show()
    app.exec()
    
    