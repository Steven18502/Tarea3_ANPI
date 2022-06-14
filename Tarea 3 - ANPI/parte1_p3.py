import sys
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from parte1_p2 import *
import sympy as sp
from sympy import *
from sympy import sympify
from sympy.sets import Interval
from sympy.calculus.util import continuous_domain

class GodWindow(QTabWidget):
    def __init__(self, parent = None):
        super(GodWindow, self).__init__(parent)
        # Dimensiones de la ventana
        self.setMaximumSize(900, 700)
        self.setMinimumSize(900, 700)
        self.resize(900, 700)

        # Fondo de la ventana principal
        self.setStyleSheet("background-color: #787272;")

        # Crear instancias de las ventanas
        self.main_tab = MainWindow()
        self.help_tab = HelpWindow()

        # Añadir las ventanas a los Tabs
        self.addTab(self.main_tab,"Tab 1")
        self.addTab(self.help_tab,"Tab 2")

        # Desplegar la configuracion de los tabs
        self.customize_mainUI()
        self.customize_helpUI()

        self.setWindowTitle("Calculadora de Integrales Definidas")
        self.tabBar().hide()


    def open_help_tab(self):
        self.setCurrentWidget(self.help_tab)

    def open_main_tab(self):
        self.setCurrentWidget(self.main_tab)

    def customize_mainUI(self):
        # Configurar el botón 1
        self.main_tab.help_button.clicked.connect(self.open_help_tab)

    def customize_helpUI(self):
        # Configurar el botón 1
        self.help_tab.main_button.clicked.connect(self.open_main_tab)

class MainWindow(QWidget):
    def __init__(self, parent = None):
        super(MainWindow, self).__init__(parent)
        # Dimensiones de la ventana
        self.setMaximumSize(900, 700)
        self.setMinimumSize(900, 700)
        self.resize(900, 700)
        # Variable del método seleccionado
        self.metodo_seleccionado = ''
        # Frame superior
        #******************************************************************************
        self.up_widget = QFrame()
        self.up_widget.setStyleSheet("border: 0px;")
        self.up_widget.setFrameShape(QFrame.StyledPanel)
        self.up_widget.setMaximumSize(800, 70)
        self.up_widget.setMinimumSize(800, 70)
        tittle_hbox = QHBoxLayout()
        tittle_label = QLabel()
        tittle_label.setText("Calculadora de Integrales Definidas")
        tittle_label.setAlignment(Qt.AlignCenter)
        tittle_label.setFont(QFont('Greek', 28))
        tittle_label.setStyleSheet("QLabel { color : #white; }")
        tittle_hbox.addWidget(tittle_label)
        self.up_widget.setLayout(tittle_hbox)

        #******************************************************************************
        # Frame del formulario superior
        self.data_widget = QFrame()
        self.data_widget.setStyleSheet("background-color: #a7a2a2;"
                                                "color: white;")
        self.data_widget.setFrameShape(QFrame.StyledPanel)
        self.data_widget.setMaximumSize(800, 150)
        self.data_widget.setMinimumSize(800, 150)
        fbox1 = QFormLayout()
        hbox0 = QHBoxLayout()
        f_label = QLabel('f(x) = ')
        f_label.setFont(QFont('Times', 12))
        self.f_entry = QLineEdit()
        self.f_entry.setAlignment(Qt.AlignCenter)
        self.f_entry.setStyleSheet("QLineEdit {"
                                   "width: 50%;"
                                   "padding: 12px 20px;"
                                   "box-sizing: border-box;"
                                   "font-size: 18px;"
                                   "background-color: #cdc3c3;"
                                   "}")

        hbox0.addWidget(f_label)
        hbox0.addWidget(self.f_entry)
        fbox1.addRow(hbox0)
        fbox1.setVerticalSpacing(50)
        hbox1 = QHBoxLayout()
        # Validador de los entrys
        a_label = QLabel('a = ')
        a_label.setFont(QFont('Times', 12))
        self.a_entry = QLineEdit()
        self.a_entry.setAlignment(Qt.AlignCenter)
        self.a_entry.setStyleSheet("font-size: 18px;"
                                   "background-color: #cdc3c3;"
                                   "width: 10px;")
        b_label = QLabel('b = ')
        b_label.setFont(QFont('Times', 12))
        self.b_entry = QLineEdit()
        self.b_entry.setAlignment(Qt.AlignCenter)
        self.b_entry.setStyleSheet("font-size: 18px;"
                                   "background-color: #cdc3c3;"
                                   "width: 10px;")
        hbox1.addWidget(a_label)
        hbox1.addWidget(self.a_entry)
        hbox1.addStretch()
        hbox1.addStretch()
        hbox1.addStretch()
        hbox1.addWidget(b_label)
        hbox1.addWidget(self.b_entry)
        fbox1.addRow(hbox1)
        self.data_widget.setLayout(fbox1)

        #******************************************************************************
        self.methods_widget = QFrame()
        self.methods_widget.setStyleSheet("background-color: #a7a2a2;"
                                          "color: white;"
                                          "border: 0px;")
        self.methods_widget.setFrameShape(QFrame.StyledPanel)
        vbox2 = QVBoxLayout()
        self.methods_widget.tabs = QTabWidget()
        self.methods_widget.tabs.setStyleSheet("background-color: #a7a2a2;"
                                                     "font-size: 14px;")

        self.methods_widget.setStyleSheet('''
        QTabWidget::tab-bar {
            alignment: center;
        }''')
        self.methods_widget.ms_tab = QFrame()
        self.methods_widget.ms_tab.setStyleSheet("color: white;")
        self.methods_widget.mc_tab = QFrame()
        self.methods_widget.mc_tab.setStyleSheet("color: white;")
        self.methods_widget.tabs.addTab(self.methods_widget.ms_tab, "Métodos Simples")
        self.methods_widget.tabs.addTab(self.methods_widget.mc_tab, "Métodos Compuestos")
        # Métodos Simples tab
        hbox2 = QVBoxLayout()
        
        self.radioGroup = QButtonGroup()
        self.t = QRadioButton("Trapecio")
        self.t.setFont(QFont('Times', 12))
        self.s = QRadioButton("Simpson")
        self.s.setFont(QFont('Times', 12))
        self.r = QRadioButton("Regla de Boole")
        self.r.setFont(QFont('Times', 12))
        self.t.toggled.connect(self.select_radio_button)
        self.s.toggled.connect(self.select_radio_button)
        self.r.toggled.connect(self.select_radio_button)
        hbox2.addWidget( self.t)
        hbox2.addWidget( self.s)
        hbox2.addWidget( self.r)
        self.methods_widget.ms_tab.setLayout(hbox2)

        # Métodos Compuestos Tab
        fbox2 = QFormLayout()
        hbox3 = QHBoxLayout()

        # Radio Group para enlazar todos los radio buttons
        self.t2 = QRadioButton("Trapecio Compuesto")
        self.t2.setFont(QFont('Times', 12))
        self.s2 = QRadioButton("Simpson Compuesto")
        self.s2.setFont(QFont('Times', 12))
        self.r2 = QRadioButton("Cuadraturas Gaussianas")
        self.r2.setFont(QFont('Times', 12))
        self.t2.toggled.connect(self.select_radio_button)
        self.s2.toggled.connect(self.select_radio_button)
        self.r2.toggled.connect(self.select_radio_button)
        hbox3.addWidget( self.t2)
        hbox3.addWidget( self.s2)
        hbox3.addWidget( self.r2)
        
        # Enlazar los radio buttons al grupo
        self.radioGroup.addButton(self.t)
        self.radioGroup.addButton(self.s)
        self.radioGroup.addButton(self.r)
        self.radioGroup.addButton(self.t2)
        self.radioGroup.addButton(self.s2)
        self.radioGroup.addButton(self.r2)
        
        fbox2.setVerticalSpacing(50)
        fbox2.addRow(hbox3)
        fbox2.setVerticalSpacing(50)
        hboxc = QHBoxLayout()
        puntos_label = QLabel('Puntos a utilizar = ')
        puntos_label.setFont(QFont('Times', 12))
        self.puntos_entry = QLineEdit()
        self.puntos_entry.setAlignment(Qt.AlignCenter)
        self.puntos_entry.setStyleSheet("font-size: 18px;"
                                        "background-color: #cdc3c3;"
                                        "width: 10px;")
        self.puntos_entry.setEnabled(False)
        hboxc.addWidget(puntos_label)
        hboxc.addWidget(self.puntos_entry)
        hboxc.addStretch()
        hboxc.addStretch()
        hboxc.addStretch()
        fbox2.addRow(hboxc)
        self.methods_widget.mc_tab.setLayout(fbox2)

        # Añadir los tabs al layout
        vbox2.addWidget(self.methods_widget.tabs)
        self.methods_widget.setLayout(vbox2)

        #******************************************************************************
        self.down_widget = QFrame()
        self.down_widget.setStyleSheet("background-color: #a7a2a2;")
        self.down_widget.setFrameShape(QFrame.StyledPanel)
        self.down_widget.setMaximumSize(800, 100)
        self.down_widget.setMinimumSize(800, 100)

        down_hbox = QHBoxLayout()
        self.help_button = QPushButton("Ayuda")
        self.help_button.setStyleSheet("QPushButton {"
                                      "background-color: #6b6767;"
                                      "font-size: 14px;"
                                      "padding: 10px 20px;"
                                      "color: white;"
                                      "border: 3px solid #000000;"
                                      "}"
                                      "QPushButton:hover {"
                                      "background-color: #c2bdbd;"
                                      "}")

        self.calcular_button = QPushButton("Calcular")
        self.calcular_button.setStyleSheet("QPushButton {"
                                           "background-color: #65806a;"
                                           "font-size: 14px;"
                                           "padding: 15px 30px;"
                                           "color: black;"
                                           "border: 3px solid #FFFFFF;"
                                           "}"
                                           "QPushButton:hover {"
                                           "background-color: #148830;"
                                           "}")

        self.calcular_button.clicked.connect(self.calcular)

        self.reset_button = QPushButton("Borrar")
        self.reset_button.setStyleSheet("QPushButton {"
                                           "background-color: #6b6767;"
                                           "font-size: 14px;"
                                           "padding: 10px 20px;"
                                           "color: white;"
                                           "border: 3px solid #000000;"
                                           "}"
                                           "QPushButton:hover {"
                                           "background-color: #c2bdbd;"
                                           "}")

        self.reset_button.clicked.connect(self.reset)
        self.resultado_label = QLabel('')
        self.resultado_label.setFont(QFont('Times', 14))
        self.resultado_label.setStyleSheet("QLabel { color : white; }")
        self.error_label = QLabel('')
        self.error_label.setFont(QFont('Times', 14))
        self.error_label.setStyleSheet("QLabel { color : white; }")
        down_hbox.addWidget(self.calcular_button)
        down_hbox.addWidget(self.resultado_label)
        down_hbox.addStretch()
        down_hbox.addWidget(self.error_label)
        down_hbox.addStretch()
        down_hbox.addWidget(self.reset_button)
        down_hbox.addWidget(self.help_button)
        self.down_widget.setLayout(down_hbox)
        self.resultado_label.setText('Aproximación: ')
        self.error_label.setText('Error: ')


        # Splitter
        splitter2 = QSplitter(Qt.Vertical)
        splitter2.setStyleSheet("border: 0px;")
        splitter2.addWidget(self.up_widget)
        splitter2.addWidget(self.data_widget)
        splitter2.addWidget(self.methods_widget)
        splitter2.addWidget(self.down_widget)


        hbox = QHBoxLayout(self)
        hbox.addWidget(splitter2)
        self.setLayout(hbox)

    def select_radio_button(self):
        radioBtn = self.sender()
        seleccionado = radioBtn.text()
        self.metodo_seleccionado = seleccionado
        if seleccionado == "Trapecio" or seleccionado == "Simpson" or seleccionado == "Regla de Boole":
            self.puntos_entry.setEnabled(False)
            self.puntos_entry.setText('')
        else:
            self.puntos_entry.setEnabled(True)


    def reset(self):
        # Limpiar todos los campos de texto
        self.f_entry.setText('')
        self.a_entry.setText('')
        self.b_entry.setText('')
        self.puntos_entry.setText('')
        self.resultado_label.setText('Aproximación:')
        self.error_label.setText('Error: ')
        # Deseleccionar todos los radio buttons
        self.radioGroup.setExclusive(False)
        self.s.setChecked(False)
        self.t.setChecked(False)
        self.r.setChecked(False)
        self.s2.setChecked(False)
        self.t2.setChecked(False)
        self.r2.setChecked(False)
        self.radioGroup.setExclusive(True)
        # Deshabilitar el entry de la cantidad de puntos
        self.puntos_entry.setEnabled(False)
        # Limpiar la variable del método seleccionado
        self.metodo_seleccionado = ''

    def validador(self, f, a, b):
        x = sp.symbols('x')
        try:
            expresion = sympify(f)
            funciones = list(expresion.atoms(Function))
            simbolos = list(expresion.free_symbols)
            dominio = continuous_domain(expresion, x, Interval(a, b))
            dominio_original = Interval(a,b)


            if(validar_simbolos(simbolos) == False):
                QMessageBox.critical(self, "Error",
                                     "La funcion está mal escrita. Puede que haya símbolos no permitidos"
                                     " o mal utilizados.")
                return False
            if(validar_funciones(funciones) == False):
                QMessageBox.critical(self, "Error",
                                     "La funcion está mal escrita. Puede que haya funciones no permitidas"
                                     " o mal utilizadas.")
                return False
            if(dominio != dominio_original):
                QMessageBox.critical(self, "Error",
                                     "La función no es contínua en el intervalo seleccionado.")
                return False

            else:
                return True
        except:
            QMessageBox.critical(self, "Error",
                                 "La función está mal escrita o tiene símbolos/funciones no disponibles actualmente.")
        return False

    def calcular(self):
        # Verificar que los primeros campos no sean vacios
        f = self.f_entry.text()
        a = self.a_entry.text()
        b = self.b_entry.text()
        puntos = self.puntos_entry.text()
        metodo = self.metodo_seleccionado
        signo = 1
        # Verificar que los primeros campos de texto tengan informacion
        if (f != '') and (a != '') and (b != ''):
            # Verificar que se haya seleccionado un método
            if (metodo != ''):
                # Validar que se haya ingresado un intervalo numérico
                try:
                    a = float(a)
                    b = float(b)
                except Exception:
                    QMessageBox.warning(self, "Advertencia", "El intervalo debe ser numérico.")
                    return

                # Validar el orden del intervalo
                if (a > b):
                    aux = a
                    a = b
                    b = aux
                    signo = -1

                # Validar si el intervalo es el mismo número
                if (a == b):
                    self.resultado_label.setText('Aproximación: ' + str(0))
                    self.error_label.setText('Error: ' + str(0))
                    return

                # Validar la estructura de la función
                if (self.validador(f, a, b) == False):
                    return

                if (metodo == 'Trapecio'):
                    s = trapecio_simple(f, a, b)
                    I = s[0]
                    er = s[1]
                    try:
                        I = float("{:.8f}".format(s[0]))
                        er = float("{:.8f}".format(s[1]))
                    except Exception:
                        pass

                    self.resultado_label.setText('Aproximación: ' + str(I*signo))
                    self.error_label.setText('Error: ' + str(er))

                elif (metodo == 'Simpson'):
                    s = simpson_simple(f, a, b)
                    I = s[0]
                    er = s[1]
                    try:
                        I = float("{:.8f}".format(s[0]))
                        er = float("{:.8f}".format(s[1]))
                    except Exception:
                        pass
                    self.resultado_label.setText('Aproximación: ' + str(I*signo))
                    self.error_label.setText('Error: ' + str(er))

                elif (metodo == 'Regla de Boole'):
                    s = boole(f, a, b)
                    I = s[0]
                    er = s[1]
                    try:
                        I = float("{:.8f}".format(s[0]))
                        er = float("{:.8f}".format(s[1]))
                    except Exception:
                        pass
                    self.resultado_label.setText('Aproximación: ' + str(I*signo))
                    self.error_label.setText('Error: ' + str(er))

                elif (metodo == 'Trapecio Compuesto'):
                    if (puntos != ''):
                        n = int(puntos)
                        s = trapecio_compuesto(f, a, b, n)
                        I = float("{:.8f}".format(s[0]))
                        er = float("{:.8f}".format(s[1]))
                        self.resultado_label.setText('Aproximación: ' + str(I*signo))
                        self.error_label.setText('Error: ' + str(er))
                    else:
                        QMessageBox.warning(self, "Advertencia", "No se ha indicado la cantidad de puntos a utilizar.")

                elif (metodo == 'Simpson Compuesto'):
                    if (puntos != ''):
                        n = int(puntos)
                        s = simpson_compuesto(f, a, b, n)
                        I = float("{:.8f}".format(s[0]))
                        er = float("{:.8f}".format(s[1]))
                        self.resultado_label.setText('Aproximación: ' + str(I*signo))
                        self.error_label.setText('Error: ' + str(er))
                    else:
                        QMessageBox.warning(self, "Advertencia", "No se ha indicado la cantidad de puntos a utilizar.")

                elif (metodo == 'Cuadraturas Gaussianas'):
                    if (puntos != ''):
                        n = int(puntos)
                        if (n < 5):
                            s = cuadraturas_gaussianas(f, a, b, n)
                            I = float("{:.8f}".format(s[0]))
                            er = float("{:.8f}".format(s[1]))
                            self.resultado_label.setText('Aproximación: ' + str(I*signo))
                            self.error_label.setText('Error: ' + str(er))
                        else:
                            QMessageBox.warning(self, "Advertencia: Cuadraturas Gaussianas", "La cantidad de puntos debe ser menor a 5 cuando"
                                                                                             " se usen las cuadraturas"
                                                                                             " gaussianas.")

                    else:
                        QMessageBox.warning(self, "Advertencia", "No se ha indicado la cantidad de puntos a utilizar.")
            else:
                QMessageBox.warning(self, "Advertencia", "No se ha seleccionado ningún método.")
        else:
                QMessageBox.warning(self, "Advertencia", "Faltan campos por llenar")


class HelpWindow(QWidget):
    def __init__(self, parent = None):
        super(HelpWindow, self).__init__(parent)
        # Dimensiones de la ventana
        self.setMaximumSize(900, 700)
        self.setMinimumSize(900, 700)
        self.resize(900, 700)

        # Frame superior
        #******************************************************************************
        self.up_widget = QFrame()
        self.up_widget.setStyleSheet("border: 0px;")
        self.up_widget.setFrameShape(QFrame.StyledPanel)
        self.up_widget.setMaximumSize(800, 70)
        self.up_widget.setMinimumSize(800, 70)

        tittle_hbox = QHBoxLayout()
        tittle_label = QLabel()
        tittle_label.setText("Ayuda")
        tittle_label.setAlignment(Qt.AlignCenter)
        tittle_label.setFont(QFont('Greek', 28))
        tittle_label.setStyleSheet("QLabel { color : #0C1B33; }")

        tittle_hbox.addWidget(tittle_label)
        self.up_widget.setLayout(tittle_hbox)

        #******************************************************************************
        # Frame central
        self.center_widget = QFrame()
        self.center_widget.setStyleSheet("background-color: #a7a2a2;"
                                                "color: white;")
        self.center_widget.setFrameShape(QFrame.StyledPanel)
        self.center_widget.setMaximumSize(800, 500)
        self.center_widget.setMinimumSize(800, 500)

        # Definiendo el layout contenedor del scroll
        vbox_layout = QVBoxLayout()
        groupBox = QGroupBox()

        # Añadiendo el contenido del Scroll Area
        l1 = QLabel()
        l1.setText("Paso 1: Ingresar la función f(x) en el primer espacio que \n"
                   "sale debajo del título de la aplicación \n"
                   "Después adregar el intervalo [a,b] en donde se desea hacer \n"
                   "la aproximación de la integral definida")
        l1.setAlignment(Qt.AlignCenter)
        l1.setFont(QFont('Times', 14))
        l1.setStyleSheet("QLabel { color : white; }")

        l2 = QLabel()
        l2.setText("Paso 2: Seleccionar un método numérico \n"
                   "Puedes elegir entre los métodos simples para resolver:  \n"
                   "   - Trapecio  \n"
                   "   - Simpson  \n"
                   "   - Regla de Boole  \n"
                   "O elegir entre los métodos Complejos: \n"
                   "   - Trapecio  Compuesto \n"
                   "   - Simpson Compuesto \n"
                   "   - Regla de Boole Compuesto \n"
                   " \n"
                   "Si eliges métodos compuestos es necesario indicar la cantidad \n"
                   "de puntos a utilizar para hacer el cálculo \n"
                   )
        l2.setAlignment(Qt.AlignCenter)
        l2.setFont(QFont('Times', 14))
        l2.setStyleSheet("QLabel { color : white; }")

        l3 = QLabel()
        l3.setText("Paso 3: Calcular la integral definida \n"
                   "Si todos los valores son ingresados correctamente \n"
                   "El sistema te mostrará el resultado approximado \n "
                   "y el error obetenido con el resultado")
        l3.setAlignment(Qt.AlignCenter)
        l3.setFont(QFont('Times', 14))
        l3.setStyleSheet("QLabel { color : white; }")

        l4 = QLabel()
        l4.setText("Puedes borrar los datos ingresados con el botón de borrar.")
        l4.setAlignment(Qt.AlignCenter)
        l4.setFont(QFont('Times', 14))
        l4.setStyleSheet("QLabel { color : white; }")

        

        espacio_label = QLabel('')
        espacio_label.setStyleSheet("margin: 30px 50px;")

        vbox_layout.addWidget(l1)
        vbox_layout.addWidget(espacio_label)
        vbox_layout.addWidget(l2)
        vbox_layout.addWidget(espacio_label)
        vbox_layout.addWidget(l3)
        vbox_layout.addWidget(espacio_label)
        vbox_layout.addWidget(l4)
        vbox_layout.addWidget(espacio_label)

        groupBox.setLayout(vbox_layout)
        scroll = QScrollArea()
        scroll.setWidget(groupBox)
        scroll.setWidgetResizable(True)
        scroll.setFixedHeight(470)
        layout = QVBoxLayout()
        layout.addWidget(scroll)
        self.center_widget.setLayout(layout)

        #******************************************************************************
        # Frame Inferior
        self.down_widget = QFrame()
        self.down_widget.setStyleSheet("background-color: #a7a2a2;"
                                         "color: white;"
                                       "border: 0px;")
        self.down_widget.setFrameShape(QFrame.StyledPanel)
        self.down_widget.setMaximumSize(800, 70)
        self.down_widget.setMinimumSize(800, 70)
        self.main_button = QPushButton("Atras")
        self.main_button.resize(60, 60)
        self.main_button.setStyleSheet("QPushButton {"
                                       "background-color: #6b6767;"
                                       "font-size: 14px;"
                                       "padding: 15px 15px;"
                                       "color: black;"
                                       "border: 3px solid black;"
                                       "}"
                                       "QPushButton:hover {"
                                       "background-color: #c2bdbd;"
                                       "}")

        down_hbox = QHBoxLayout()
        down_hbox.addStretch(10)
        down_hbox.addWidget(self.main_button)
        down_hbox.addStretch(10)
        self.down_widget.setLayout(down_hbox)

        # Configuracion de los splitters
        splitter = QSplitter(Qt.Vertical)
        splitter.addWidget(self.up_widget)
        splitter.addWidget(self.center_widget)
        splitter.addWidget(self.down_widget)


        main_hbox = QHBoxLayout(self)
        main_hbox.addWidget(splitter)
        self.setLayout(main_hbox)


def main():
    app = QApplication(sys.argv)
    ex = GodWindow()
    ex.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
