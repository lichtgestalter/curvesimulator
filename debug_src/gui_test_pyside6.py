import sys
from PySide6.QtWidgets import QApplication, QMainWindow, QPushButton, QLabel, QVBoxLayout, QWidget
from PySide6.QtCore import Qt

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("PySide6 Test")
        self.setFixedSize(300, 200)

        # Central widget und Layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)

        # Label
        self.label = QLabel("Hallo PySide6!")
        layout.addWidget(self.label)

        # Button mit Event
        button = QPushButton("Klick mich!")
        button.clicked.connect(self.on_button_click)
        layout.addWidget(button)

    def on_button_click(self):
        self.label.setText("Button geklickt! 😊")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
