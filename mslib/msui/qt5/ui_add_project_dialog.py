# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ui_add_project.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_addProjectDialog(object):
    def setupUi(self, addProjectDialog):
        addProjectDialog.setObjectName("addProjectDialog")
        addProjectDialog.resize(437, 212)
        self.buttonBox = QtWidgets.QDialogButtonBox(addProjectDialog)
        self.buttonBox.setGeometry(QtCore.QRect(60, 150, 301, 32))
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.label = QtWidgets.QLabel(addProjectDialog)
        self.label.setGeometry(QtCore.QRect(90, 50, 59, 15))
        self.label.setObjectName("label")
        self.label_2 = QtWidgets.QLabel(addProjectDialog)
        self.label_2.setGeometry(QtCore.QRect(70, 80, 71, 20))
        self.label_2.setObjectName("label_2")
        self.path = QtWidgets.QLineEdit(addProjectDialog)
        self.path.setGeometry(QtCore.QRect(150, 50, 211, 23))
        self.path.setObjectName("path")
        self.description = QtWidgets.QTextEdit(addProjectDialog)
        self.description.setGeometry(QtCore.QRect(150, 80, 211, 41))
        self.description.setObjectName("description")

        self.retranslateUi(addProjectDialog)
        self.buttonBox.accepted.connect(addProjectDialog.accept)
        self.buttonBox.rejected.connect(addProjectDialog.reject)
        QtCore.QMetaObject.connectSlotsByName(addProjectDialog)

    def retranslateUi(self, addProjectDialog):
        _translate = QtCore.QCoreApplication.translate
        addProjectDialog.setWindowTitle(_translate("addProjectDialog", "Add project"))
        self.label.setText(_translate("addProjectDialog", "path"))
        self.label_2.setText(_translate("addProjectDialog", "description"))

