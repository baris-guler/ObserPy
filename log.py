from PyQt5.QtWidgets import QDialog, QVBoxLayout, QLabel, QLineEdit, QPushButton, QTableWidgetItem, QFileDialog, QMessageBox, QTableWidget, QHeaderView
from PyQt5.QtCore import Qt, QRunnable, QThreadPool, pyqtSignal, QObject, pyqtSlot, QDateTime, QDate, QTime
from PyQt5 import QtWidgets, uic
from PyQt5.QtCore import QDate, Qt
import os
import json
import shutil
import sys
from PyQt5.QtWidgets import QApplication
import subprocess
from PIL import Image

log = None
logedit = None
selected_date = None
mapf = None
lcf = None
edit_mode = False

def open_photo(file_path):
    if file_path.endswith((".jpg", ".jpeg", ".png")):
        subprocess.run(["xdg-open", file_path])
    else:
        print("Unsupported file format.")

def get_subdir_names(directory):
    dir_names = []
    for root, dirs, files in os.walk(directory):
        for dir in dirs:
            dir_names.append(dir)
    return dir_names

def open_log(star_name):
    global log
    log = uic.loadUi("log.ui")
    log_note_path = f'log/{star_name}/{star_name}_note.txt'
    if os.path.exists(log_note_path):
        log.log_note.setText(open_note(log_note_path))
    log.log_note.setDisabled(True)
    log.log_save.clicked.connect(lambda: save_note(star_name))
    log.log_save.setDisabled(True)
    log.log_cancel.setDisabled(True)
    log.log_cancel.clicked.connect(lambda: cancel_note(star_name))
    log_label = log.log_label
    log_label.setText(star_name)
    log_add = log.log_add
    log_add.clicked.connect(lambda: add_log(star_name))
    log_edit = log.log_edit
    log_edit.setDisabled(True)
    log_del = log.log_del
    log_del.setDisabled(True)
    header = log.log_table.horizontalHeader()
    header.setSectionResizeMode(0, QHeaderView.Stretch)
    log.log_table.cellClicked.connect(date_selected)
    log.log_table.setColumnCount(2)
    log.log_table.itemDoubleClicked.connect(lambda: display_log(None, None))
    log.log_editbox.stateChanged.connect(editbox_changed)
    log.log_display.clicked.connect(lambda: display_log(None, None))
    
    header_labels = ["Date", "Telescope"]
    log.log_table.setHorizontalHeaderLabels(header_labels)

    directory_path = f'log/{star_name}'
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
    subdir_names = get_subdir_names(directory_path)
    
    log.log_table.setRowCount(len(subdir_names))  # Set the number of rows
    
    # Populate the table with file names
    for row, folder_name in enumerate(subdir_names):
        item = QTableWidgetItem(folder_name)
        item.setTextAlignment(Qt.AlignCenter) 
        log.log_table.setItem(row, 0, item)
        
        #open folder and open the first file that ends with .json
        folder_path = f'log/{star_name}/{folder_name}'
        json_files = [pos_json for pos_json in os.listdir(folder_path) if pos_json.endswith('.json')]
        if len(json_files) > 0:
            json_file = json_files[0]
            json_path = f'log/{star_name}/{folder_name}/{json_file}'
            with open(json_path, 'r') as f:
                json_data = json.load(f)
                log.log_table.setItem(row, 1, QTableWidgetItem(json_data['telescope']))

    log.show()
    log.activateWindow()
    
def display_log(row, column):
    star_name = log.log_label.text()
    if row is None:
        row = log.log_table.currentRow()
    log_date = log.log_table.item(row,0).text()
    edit_log(star_name, log_date, display=True)
    
def open_note(file_path):
    with open(file_path, 'r') as f:
        note = f.read()
    return note
    
def save_note(star_name):
    note = log.log_note.toPlainText()
    file_path = f'log/{star_name}/{star_name}_note.txt'
    with open(file_path, 'w') as f:
        f.write(note)
    log.log_editbox.setChecked(False)
    
def cancel_note(star_name):
    log_note_path = f'log/{star_name}/{star_name}_note.txt'
    log.log_note.setText(open_note(log_note_path))
    log.log_editbox.setChecked(False)
    
def editbox_changed(state):
    if state == 2:
        log.log_note.setDisabled(False)
        log.log_save.setDisabled(False)
        log.log_cancel.setDisabled(False)
    else:
        log.log_note.setDisabled(True)
        log.log_save.setDisabled(True)
        log.log_cancel.setDisabled(True)
    
def date_selected(row, column):
    global selected_date
    log.log_table.selectRow(row)
    item = log.log_table.item(row, 0)
    if item is not None:
        text = item.text()
        selected_date = text
        star_name = log.log_label.text()
        log.log_edit.setDisabled(False)
        log.log_edit.clicked.connect(lambda: edit_log(star_name, selected_date))
        log.log_del.setDisabled(False)
        log.log_del.clicked.connect(lambda: del_log(star_name, selected_date))
        
def add_log(star_name):
    global logedit
    global edit_mode
    edit_mode = False
    logedit = uic.loadUi("log_edit.ui")
    logedit_star_name = logedit.logedit_star_name
    logedit_star_name.setText(star_name)
    logedit_star_name.setDisabled(True)
    logedit_filter = logedit.logedit_filter
    logedit_exp_time = logedit.logedit_exp_time
    logedit_note = logedit.logedit_note
    logedit_binning = logedit.logedit_binning
    logedit_lc_edit = logedit.logedit_lc_edit
    logedit_lc_edit.clicked.connect(lambda: lc_edit(star_name))
    logedit_lc_display = logedit.logedit_lc_display
    logedit_lc_display.setDisabled(True)
    logedit_lc_display.clicked.connect(lambda: lc_display())
    logedit_map_edit = logedit.logedit_map_edit
    logedit_map_edit.clicked.connect(lambda: map_edit(star_name))
    logedit_map_display = logedit.logedit_map_display
    logedit_map_display.setDisabled(True)
    logedit_map_display.clicked.connect(lambda: map_display())
    logedit_date = logedit.logedit_date
    logedit_date.setDisabled(False)
    current_date = QDate.currentDate()
    yesterday = current_date.addDays(-1)
    logedit_date.setDate(yesterday)
    logedit_save = logedit.logedit_save
    logedit_save.clicked.connect(lambda: save_log(star_name))
    logedit_cancel = logedit.logedit_cancel
    logedit_cancel.clicked.connect(lambda: cancel_log(star_name))
    logedit.show()
    logedit.activateWindow()

def edit_log(star_name, selected_date, display=False):
    global logedit
    global map_path
    global lc_path
    global edit_mode
    
    edit_mode = True
    logedit = uic.loadUi("log_edit.ui")
    file_path = "log/" + star_name + "/" + selected_date + "/" + selected_date + "_log.json"
    map_path = "log/" + star_name + "/" + selected_date + "/" + selected_date + "_map.png"
    lc_path = "log/" + star_name + "/" + selected_date + "/" + selected_date + "_lc.png"
    
    if display:
        logedit.logedit_save.setDisabled(True)
    else:
        logedit.logedit_save.setDisabled(False)
    
    with open(file_path, "r") as file:
        json_data = file.read()
        log_dict = json.loads(json_data)

    # Now you can access the dictionary values as needed
    star_name = log_dict["star_name"]
    filter_value = log_dict["filter"]
    exp_time = log_dict["exp_time"]
    note = log_dict["note"]
    date = log_dict["date"]
    binning = log_dict["binning"]
    telescope = log_dict["telescope"]
    logedit.logedit_star_name.setText(star_name)
    logedit.logedit_star_name.setDisabled(True)
    logedit.logedit_filter.setText(filter_value)
    logedit.logedit_exp_time.setText(exp_time)
    logedit.logedit_note.setText(note)
    logedit.logedit_date.setDate(QDate.fromString(date, "yyyy-MM-dd"))
    logedit.logedit_date.setDisabled(True)
    logedit.logedit_binning.setText(binning)
    logedit.logedit_telescope.setText(telescope)
    logedit.logedit_save.clicked.connect(lambda: save_log(star_name))
    logedit.logedit_lc_edit.clicked.connect(lambda: lc_edit(star_name))
    logedit_lc_display = logedit.logedit_lc_display
    if lcf == None and os.path.exists(lc_path) == False:
        logedit_lc_display.setDisabled(True)
    else:
        logedit_lc_display.setDisabled(False)
    logedit_lc_display.clicked.connect(lambda: lc_display())
    logedit_map_edit = logedit.logedit_map_edit
    logedit_map_edit.clicked.connect(lambda: map_edit(star_name))
    logedit_map_display = logedit.logedit_map_display
    if mapf == None and os.path.exists(map_path) == False:
        logedit_map_display.setDisabled(True)
    else:
        logedit_map_display.setDisabled(False)
    logedit_map_display.clicked.connect(lambda: map_display())
    logedit_cancel = logedit.logedit_cancel
    logedit_cancel.clicked.connect(lambda: cancel_log(star_name))
    logedit.show()
    logedit.activateWindow()

def del_log(star_name, selected_date):
    confirm_deletion(star_name, selected_date)

def confirm_deletion(star_name, selected_date):
    
    message_box = QMessageBox()
    message_box.setIcon(QMessageBox.Warning)
    message_box.setWindowTitle("Confirmation")
    message_box.setText("Are you sure you want to delete the log?")
    message_box.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
    message_box.setDefaultButton(QMessageBox.No)

    reply = message_box.exec_()
    if reply == QMessageBox.Yes:
        folder_path = "log/" + star_name + "/" + selected_date
        if folder_path:
            shutil.rmtree(folder_path)
            open_log(star_name)
    if reply == QMessageBox.No:
        message_box.close()
        
def lc_edit(star_name):
    global lcf
    options = QFileDialog.Options()
    options |= QFileDialog.DontUseNativeDialog
    
    file_dialog = QFileDialog()
    file_dialog.setOptions(options)
    
    file_dialog.setFileMode(QFileDialog.ExistingFile)
    
    if file_dialog.exec_():
        selected_files = file_dialog.selectedFiles()
        lc_file = selected_files[0]
        file_extension = os.path.splitext(lc_file)[1].lower()
        if file_extension not in ['.png', '.jpg', '.jpeg']:
            message_box = QMessageBox()
            message_box.setIcon(QMessageBox.Critical)
            message_box.setWindowTitle("Error")
            message_box.setText("File must be .png or .jpg file")
            message_box.exec_()
        else:
            lcf = lc_file
            logedit.logedit_lc_display.setDisabled(False)

def lc_display():
    global lc_path
    if lcf != None:
        im = Image.open(lcf)
        im.show()
    elif lc_path.endswith((".jpg", ".jpeg", ".png")):
        im = Image.open(lc_path)
        im.show()
    else:
        print("Unsupported file format.")  

def map_edit(star_name):
    global mapf
    options = QFileDialog.Options()
    options |= QFileDialog.DontUseNativeDialog 
    
    file_dialog = QFileDialog()
    file_dialog.setOptions(options)
    
    file_dialog.setFileMode(QFileDialog.ExistingFile)
    
    if file_dialog.exec_():
        selected_files = file_dialog.selectedFiles()
        map_file = selected_files[0]
        file_extension = os.path.splitext(map_file)[1].lower()
        if file_extension not in ['.png', '.jpg', '.jpeg']:
            message_box = QMessageBox()
            message_box.setIcon(QMessageBox.Critical)
            message_box.setWindowTitle("Error")
            message_box.setText("File must be .png or .jpg file")
            message_box.exec_()
        else:
            mapf = map_file
            logedit.logedit_map_display.setDisabled(False)

def map_display():
    global map_path
    if mapf != None:
        im = Image.open(mapf)
        im.show()
    elif map_path.endswith((".jpg", ".jpeg", ".png")):
        im = Image.open(map_path)
        im.show()
    else:
        print("Unsupported file format.")   

def save_log(star_name):
    global mapf
    global lcf
    global map_path
    global lc_path
    log_dict = {}
    log_dict["star_name"] = star_name
    log_dict["filter"] = logedit.logedit_filter.text()
    log_dict["exp_time"] = logedit.logedit_exp_time.text()
    log_dict["note"] = logedit.logedit_note.toPlainText()
    log_dict["date"] = logedit.logedit_date.date().toString("yyyy-MM-dd")
    log_dict["binning"] = logedit.logedit_binning.text()
    log_dict["telescope"] = logedit.logedit_telescope.text()
    json_data = json.dumps(log_dict, indent=4)
    
    if not edit_mode:
        file_path = f"log/{star_name}/{logedit.logedit_date.date().toString('yyyy-MM-dd')}"
        if os.path.exists(file_path):
            message_box = QMessageBox()
            message_box.setIcon(QMessageBox.Warning)
            message_box.setWindowTitle("Warning")
            message_box.setText("There is already a log for this day. Do you want to replace it?")
            replace_button = QPushButton("Replace")
            create_new_button = QPushButton("Create New")
            cancel_button = QPushButton("Cancel")

            message_box.addButton(replace_button, QMessageBox.YesRole)
            message_box.addButton(create_new_button, QMessageBox.NoRole)
            message_box.addButton(cancel_button, QMessageBox.RejectRole)
            reply = message_box.exec_()
            if reply == 0:
                shutil.rmtree(file_path)
                os.mkdir(file_path)
            elif reply == 1:
                file_path_old = file_path
                i = 1
                while os.path.exists(file_path_old + "_" + str(i)):
                    i += 1
                print("sea")
                file_path = file_path_old + "_" + str(i)
                os.mkdir(file_path)
            first_selected_row_text = file_path.split("/")[-1]
        else:
            os.mkdir(file_path)
            first_selected_row_text = file_path.split("/")[-1]    
    else:
        selected_items = log.log_table.selectedItems()
        first_selected_row_text = selected_items[0].tableWidget().item(selected_items[0].row(), 0).text()
        file_path = f"log/{star_name}/{first_selected_row_text}"
        if os.path.exists(file_path):
            shutil.rmtree(file_path)
        os.mkdir(file_path)
    
    map_path = os.path.join(file_path + "/" + first_selected_row_text + "_map.png")
    if mapf is not None:
        shutil.copy2(mapf, map_path)
        mapf = None
    lc_path = os.path.join(file_path + "/" + first_selected_row_text + "_lc.png")
    if lcf is not None:
        shutil.copy2(lcf, lc_path)
        lcf = None
    
    with open(file_path + "/" + first_selected_row_text + "_log.json", "w") as file:
        file.write(json_data)
    logedit.close()
    open_log(star_name)

def cancel_log(star_name):
    global mapf
    global lcf
    logedit.close()
    mapf = None
    lcf = None
    open_log(star_name)
    
if __name__ == "__main__":
    subprocess.run(["python", "main.py"])
    