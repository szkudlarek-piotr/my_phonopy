import os

primitive_cells = []
vibrations = []
number_of_atoms = 0
supercells = []
displaced_supercells = []

class Cell:
    def __init__(self, xdim, ydim, zdim):
        self.xdim = xdim
        self.ydim = ydim
        self.zdim = zdim
        self.dict_of_atoms = {}
        print("Tworzę komorkę {} A x {} A x {} A".format(self.xdim, self.ydim, self.zdim))
    def add_atom(self, atom_obj):
        number_of_atoms = len(self.dict_of_atoms)
        new_atom_index = number_of_atoms + 1
        atom_name = "atom_" + str(new_atom_index)
        self.dict_of_atoms[atom_name] = atom_obj
    def describe_own_dimensions(self):
        print("Wymiary: {} A x {} A x {} A".format(self.xdim, self.ydim, self.zdim))
        for atom_name in self.dict_of_atoms:
            atom = self.dict_of_atoms[atom_name]
            print("{}\t{}\t{}\t{}".format(atom.element, atom.x_ang, atom.y_ang, atom.z_ang))
    def create_supercell(self, a, b, c):
        new_x_dim = self.xdim * a
        new_y_dim = self.ydim * b
        new_z_dim = self.zdim * c
        new_supercell = Cell(new_x_dim, new_y_dim, new_z_dim)
        for i in range(0, a):
            for j in range(0, b):
                for k in range(0, c):
                    starting_x = i * self.xdim
                    starting_y = j * self.ydim
                    starting_z = k * self.zdim
                    for atom_name in self.dict_of_atoms:
                        atom_obj = self.dict_of_atoms[atom_name]
                        el = atom_obj.element
                        old_x = atom_obj.x_ang
                        old_y = atom_obj.y_ang
                        old_z = atom_obj.z_ang
                        new_x = old_x + starting_x
                        new_y = old_y + starting_y
                        new_z = old_z + starting_z
                        atom_to_add = Atom(el, new_x, new_y, new_z)
                        new_supercell.add_atom(atom_to_add)
        supercells.append(new_supercell)
    def create_displaced_supercell(self, a, b, c, arr_of_displacements):
        new_x_dim = self.xdim * a
        new_y_dim = self.ydim * b
        new_z_dim = self.zdim * c
        displaced_supercell = Cell(new_x_dim, new_y_dim, new_z_dim)
        for i in range(0, a):
            for j in range(0, b):
                for k in range(0, c):
                    starting_x = i * self.xdim
                    starting_y = j * self.ydim
                    starting_z = k * self.zdim
                    sum_of_subcell_indices = i + j + k
                    if sum_of_subcell_indices % 2 == 0:
                        sign_of_displacement = 1
                    else:
                        sign_of_displacement = -1
                    #potrzebuję tego, żeby po indeksie w tablicy ustalić które drganie wziąć z tablicy
                    number_of_atom = 0
                    for atom_name in self.dict_of_atoms:
                        atom_obj = self.dict_of_atoms[atom_name]
                        displacement_for_atom = arr_of_displacements[number_of_atom]
                        x_move = displacement_for_atom[0] * sign_of_displacement
                        y_move = displacement_for_atom[1] * sign_of_displacement
                        z_move = displacement_for_atom[2] * sign_of_displacement
                        el = atom_obj.element
                        old_x = atom_obj.x_ang
                        old_y = atom_obj.y_ang
                        old_z = atom_obj.z_ang
                        new_x = old_x + starting_x
                        new_y = old_y + starting_y
                        new_z = old_z + starting_z
                        atom_to_add = Atom(el, new_x, new_y, new_z)
                        atom_to_add.displace_atom(x_move, y_move, z_move)
                        displaced_supercell.add_atom(atom_to_add)
                        number_of_atom += 1
        displaced_supercell.normalize_coordinates()
        displaced_supercell.describe_in_qe_format()
        displaced_supercells.append(displaced_supercell)
    def describe_in_qe_format(self):
        print("CELL_PARAMETERS (angstrom)")
        print("{}\t0.00000\t0.00000".format(self.xdim))
        print("0.00000\t{}\t0.00000".format(self.ydim))
        print("0.00000\t0.00000\t{}\n".format(self.zdim))
        print("ATOMIC_POSITIONS (crystal)")
        for atom_name in self.dict_of_atoms:
            atom_obj = self.dict_of_atoms[atom_name]
            x_frac = atom_obj.x_ang / self.xdim
            y_frac = atom_obj.y_ang / self.ydim
            z_frac = atom_obj.z_ang / self.zdim
            element = atom_obj.element
            print("{}\t{}\t{}\t{}".format(element, x_frac, y_frac, z_frac))
    def normalize_coordinates(self):
        for atom_name in self.dict_of_atoms:
            atom_obj = self.dict_of_atoms[atom_name]
            if atom_obj.x_ang < 0:
                atom_obj.x_ang += self.xdim
            if atom_obj.y_ang < 0:
                atom_obj.y_ang += self.ydim
            if atom_obj.z_ang < 0:
                atom_obj.z_ang += self.zdim

class Atom:
    def __init__(self, element,x_ang, y_ang, z_ang):
        self.element = element
        self.x_ang = x_ang
        self.y_ang = y_ang
        self.z_ang = z_ang
    def displace_atom(self, xvec, yvec, zvec):
        displacement_factor = 0.1
        xdisp = xvec * displacement_factor
        ydisp = yvec * displacement_factor
        zdisp = zvec * displacement_factor
        self.x_ang += xdisp
        self.y_ang += ydisp
        self.z_ang += zdisp

class QE_input:
    def __init__(self, input_path):
        self.input_path = input_path
    def read_cell(self):
        global number_of_atoms
        r = open(self.input_path, "r")
        lines = r.readlines()
        for i in range(0, len(lines)):
            line = lines[i]
            if "CELL_PARAMETERS (angstrom)" in line:
                cell_start_index = i + 1
            if "ATOMIC_POSITIONS (crystal)" in line:
                atoms_start_index = i + 1
        #rozważam tylko komórki prostopadłościenne
        x_dim = float(lines[cell_start_index].split()[0])
        y_dim = float(lines[cell_start_index+1].split()[1])
        z_dim = float(lines[cell_start_index+2].split()[2])
        my_cell = Cell(x_dim, y_dim, z_dim)
        for i in range(atoms_start_index, len(lines)):
            atom_line = lines[i]
            if len(atom_line.split()) == 4:
                element = atom_line.split()[0]
                xcoor = float(atom_line.split()[1]) * x_dim
                ycoor = float(atom_line.split()[2]) * y_dim
                zcoor = float(atom_line.split()[3]) * z_dim
                atom_object = Atom(element, xcoor, ycoor, zcoor)
                my_cell.add_atom(atom_object)
                number_of_atoms += 1
        primitive_cells.append(my_cell)

class Vibration:
    #xyz opisuje tutaj wektory nie wychylenia, a położenie drgania w komórce
    def __init__(self, freq_in_cm, x, y, z, number_in_point):
        self.freq_in_cm = freq_in_cm
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.number_in_point = number_in_point
    def add_vib_array(self, vibration_array):
        self.vib_array = vibration_array
    def tell_me_where(self):
        print("Jestem w {}_{}_{}.".format(self.x, self.y, self.z))
    def get_supercell_dims(self):
        to_show = []
        if self.x != 0:
            xdim = 2
        else:
            xdim = 1
        if self.y != 0:
            ydim = 2
        else:
            ydim = 1
        zdim = 1
        to_show.append(xdim)
        to_show.append(ydim)
        to_show.append(zdim)
        print(type(to_show))
        return to_show

class Modes_output:
    def __init__(self, path):
        self.modes_dict = {}
        r = open(path, "r")
        self.modes_lines = r.readlines()
        for i in range(2, len(self.modes_lines)):
            line = self.modes_lines[i]
            if "q =" in line:
                array_to_add = []
                x_vec = line.split()[2]
                y_vec = line.split()[3]
                z_vec = line.split()[4]
                name_of_dict = "vibrations_{}_{}_{}".format(x_vec, y_vec, z_vec)
                min_freq_in_point = float(self.modes_lines[i+2].split()[-2])
                new_vib = Vibration(min_freq_in_point, x_vec, y_vec, z_vec, 1)
                for j in range(i + 3, i + 3 + number_of_atoms):
                    vib_line = self.modes_lines[j]
                    x_vib = float(vib_line.split()[1])
                    y_vib = float(vib_line.split()[3])
                    z_vib = float(vib_line.split()[5])
                    array_to_add.append([x_vib, y_vib, z_vib])
                new_vib.add_vib_array(array_to_add)
                vibrations.append(new_vib)
    def show_modes(self):
        print(vibrations)
    def find_by_point(self):
        letters_to_point = {"G": [0, 0, 0],"GX": [0, 0.25, 0],"XG": [0, 0.25, 0], "X": [0, -0.5, 0], "XM": [0.25, -0.5,0],
        "MX": [0.25, -0.5, 0], "M": [-0.5, -0.5, 0],"MG": [0.25, 0.25, 0], "GM": [0.25, 0.25, 0], "Z": [0, 0, -0.5], "R": [0, -0.5,  -0.5],
        "A": [-0.5, -0.5,-0.5], "ZR": [0, 0.25, -0.5], "RZ": [0, 0.25, -0.5], "RA": [0.25, -0.5, -0.5], "AR": [0.25,-0.5, -0.5],
        "AZ": [0.25, 0.25, -0.5], "ZA": [0.25, 0.25, -0.5], "ZX": [0, 0.25, 0.25], "XZ": [0, 0.25, 0.25]}
        while True:
            examined_point = input("Podaj oznaczenie punktu, dla któego zwiualizować najniższe wychylenie.\n")
            if examined_point in letters_to_point:
                return letters_to_point[examined_point]
                break
            else:
                print("Nie znam tego uznaczenia, podaj inny punkt!!!")

#mother = input("Podaj ścieżkę do folderu, gdzie masz scf.in oraz matdyn.modes:\n")

mother = r"C:\Users\Piotr Szkudlarek\Desktop\doktorat\kanapki\wyniki\LF_LH_1lay\disp_35"
input_dir = os.path.join(mother, "geom.in")
modes_dir = os.path.join(mother, "matdyn.modes")
my_input = QE_input(input_dir)
my_input.read_cell()
primitive_cell = primitive_cells[0]
modes = Modes_output(modes_dir)
#modes.show_modes()
minimal_freq = 0
for i in vibrations:
    freq_value = i.freq_in_cm
    if freq_value < minimal_freq:
        minimal_freq = freq_value
for i in vibrations:
    if i.freq_in_cm == minimal_freq:
        analyzed_vib = i
        minimal_vibration = i.vib_array
        i.tell_me_where()

def get_supercell_dims(x,y,z):
    if x != 0:
        xdim = 2
    else:
        xdim = 1
    if y != 0:
        ydim = 2
    else:
        ydim = 1
    dims_array = [xdim, ydim, z]
    return dims_array

primitive_cell.create_displaced_supercell(1, 2, 1, minimal_vibration)
displaced_cell = displaced_supercells[0]
displaced_cell.describe_in_qe_format()
def displaced_from_point():
    global primitive_cell, vibrations, modes
    ar_of_vec = modes.find_by_point()
    xvec = ar_of_vec[0]
    yvec = ar_of_vec[1]
    zvec = ar_of_vec[2]
    supercell_dims = get_supercell_dims(xvec, yvec, zvec)
    x_size = supercell_dims[0]
    y_size = supercell_dims[1]

    for vib in vibrations:
        x_component = vib.x
        y_component = vib.y
        z_component = vib.z
        if x_component == xvec and y_component == yvec and z_component == z_component:
            array_to_use = vib.vib_array
            break
    print(array_to_use)
    primitive_cell.create_displaced_supercell(x_size, y_size, 1, array_to_use)

displaced_from_point()