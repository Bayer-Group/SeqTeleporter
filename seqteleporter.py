import sys
import os.path
from PyQt5.QtWidgets import QApplication, QWidget, QPushButton, QFileDialog, QVBoxLayout, QLabel

SCRIPT_DIR = os.path.dirname(os.path.abspath('__file__'))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from seqteleporter.main_libs import (generate_and_optimize_ready_to_click_modules,
                                     assemble_modules_and_generate_robot_instruction)

print(f'Here is: {SCRIPT_DIR}')


class SeqTeleporter(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.setWindowTitle('SeqTeleporter')

        self.layout = QVBoxLayout()

        self.label = QLabel('Choose an Excel file to process:')
        self.layout.addWidget(self.label)

        self.btnLoad = QPushButton('Load Excel File')
        self.btnLoad.clicked.connect(self.load_file)
        self.layout.addWidget(self.btnLoad)

        self.setLayout(self.layout)

    def load_file(self):
        options = QFileDialog.Options()
        file_name, _ = QFileDialog.getOpenFileName(
            self,
            "Select Excel File",
            "",
            "Excel Files (*.xls *.xlsx)",
            options=options
        )
        if file_name:
            self.run(file_name)

    def run(self, file_name):
        ready_to_click_modules_paths = generate_and_optimize_ready_to_click_modules(file_name)
        robot_instruction_path = assemble_modules_and_generate_robot_instruction(
            file_name,
            ready_to_click_modules_paths[0]
        )

        self.label.setText(f'Processed: {file_name}\n'
                           f'Generated: {ready_to_click_modules_paths}, \n {robot_instruction_path}')


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = SeqTeleporter()
    ex.resize(400, 200)
    ex.show()
    sys.exit(app.exec_())


