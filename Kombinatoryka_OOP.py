import itertools
import sys
import shutil
import numpy as np
import math
import string
from send2trash import send2trash


# noinspection PyShadowingNames
class Inicjalizacja(object):
    def __init__(self, mode, variance, size):
        self.mode = mode
        self.variance = variance
        self.size = size
        self.combination = None
        self.variation = None

    def variations_with_repeats(self):
        self.variation = self.variance ** self.size
        return self.variation

    def combinations_with_repeats(self):
        self.combination = int(math.factorial(self.variance + self.size - 1) / (
                math.factorial(self.size) * math.factorial(self.variance - 1)))
        return self.combination

    def __str__(self):
        if self.mode not in [2, 3, 4, 5]:
            return "Liczba możliwych polimerów: %s\nLiczba unikalnych składów: %s" % (
                self.variation, self.combination)
        elif self.mode == 2:
            return "Akcja zakończona pomyślnie"
        elif self.mode == 3:
            return "Liczba możliwych polimerów: %s" % self.variation
        elif self.mode == 4:
            return "Liczba unikalnych składów: %s" % self.combination


# noinspection PyShadowingNames
class Polymer:
    def __init__(self, sequence, variance):
        # Polimer to będzie tupple, niezmienna lista z nawiasami okrągłymi
        self.sequence = sequence
        self.variance = variance
        self.string = '-'.join(str(self.sequence).strip('()').split(', '))
        self.component_list = []
        self.zero_array = None

    def create_component_list(self):
        self.zero_array = np.zeros(self.variance)
        self.component_list = list(self.zero_array.astype(int))
        for each_monomer in self.sequence:
            self.component_list[each_monomer - 1] += 1
        return self.component_list

    def __str__(self):
        return f"Sekwencja tego polimeru to: {self.string}\nKompozycja tego polimeru to: {self.component_list}"

    def __len__(self):
        return len(self.sequence)


# noinspection PyShadowingNames
class NmerList(Inicjalizacja):
    def __init__(self, mode, variance, size):
        super().__init__(mode, variance, size)
        self.range_of_monomers = range(1, self.variance + 1)
        self.nmer_list = []

    def create_nmer_list(self):
        for xs in itertools.product(self.range_of_monomers, repeat=self.size):
            polymer = Polymer(xs, self.variance)
            self.nmer_list.append(polymer)
        return self.nmer_list

    def __len__(self):
        return len(self.nmer_list)

    def __str__(self):
        return f"Lista polimerów to: {self.nmer_list}"


# noinspection PyShadowingNames
class Matrix:
    def __init__(self, nmer_list):
        self.nmer_list = nmer_list
        self.matrix = []

    def create_visual_matrix(self):
        self.matrix = {}
        for polymer in self.nmer_list:
            polymer.create_component_list()
            if self.matrix == {}:
                self.matrix[f'{polymer.component_list}'] = [1, [polymer.string]]
            elif f'{polymer.component_list}' in self.matrix.keys():
                self.matrix[f'{polymer.component_list}'][0] += 1
                self.matrix[f'{polymer.component_list}'][1].append(polymer.string)
            else:
                self.matrix[f'{polymer.component_list}'] = [1, [polymer.string]]
        self.matrix = self.matrix.items()
        return self.matrix

    def create_visual_matrix_values(self):
        self.matrix = {}
        for polymer in self.nmer_list:
            polymer.create_component_list()
            if self.matrix == {}:
                self.matrix[f'{polymer.component_list}'] = [polymer.string]
            elif f'{polymer.component_list}' in self.matrix.keys():
                self.matrix[f'{polymer.component_list}'].append(polymer.string)
            else:
                self.matrix[f'{polymer.component_list}'] = [polymer.string]
        self.matrix = self.matrix.values()
        return self.matrix

    def create_order_matrix(self):
        for polymer in self.nmer_list:
            self.matrix.append(list(polymer.sequence))
        return self.matrix


# noinspection PyShadowingNames
class File:
    def __init__(self, variance, size, mode, matrix):
        self.variance = variance
        self.size = size
        self.mode = mode
        self.file_name = ''
        self.matrix = matrix
        self.file = None

    def map_loop_vis(self):
        a = ''
        list(map(lambda x: self.file.write(str(a + str(x) + '\n')), self.matrix))
        self.file.close()

    def map_loop_ord(self):
        self.file.write(str(self.matrix))
        self.file.close()

    def create_file_vis(self):
        self.file = open(f'{str(self.size)}-mer_{str(self.variance)}_components.txt', 'w+')
        self.file_name = f'{str(self.size)}-mer_{str(self.variance)}_components.txt'

    def create_file_ord(self):
        self.file = open(f'{str(self.size)}-mer_{str(self.variance)}_order_matrix.txt', 'w+')
        self.file_name = f'{str(self.size)}-mer_{str(self.variance)}_order_matrix.txt'

    def move_file_vis(self):
        try:
            shutil.move(self.file_name, 'E:\\Bioinformatyka\\PORT\\N-mery\\visual_info')
        except OSError:
            send2trash('E:\\Bioinformatyka\\PORT\\N-mery\\visual_info\\' + self.file_name)
            shutil.move(self.file_name, 'E:\\Bioinformatyka\\PORT\\N-mery\\visual_info')

    def move_file_ord(self):
        try:
            shutil.move(self.file_name, 'E:\\Bioinformatyka\\PORT\\N-mery\\order_matrix')
        except OSError:
            send2trash('E:\\Bioinformatyka\\PORT\\N-mery\\order_matrix\\' + self.file_name)
            shutil.move(self.file_name, 'E:\\Bioinformatyka\\PORT\\N-mery\\order_matrix')


program_start = True
positive_lower = ['y', 'yes', 'yea', 'yeah', 't', 'ta', 'tak']
positive_cap = [each.capitalize() for each in positive_lower]
positive_upper = [each.upper() for each in positive_lower]
all_p = positive_upper + positive_cap + positive_lower
negative_lower = ['n', 'no', 'nah', 'nope', 'nie', 'ni']
negative_cap = [each.capitalize() for each in negative_lower]
negative_upper = [each.upper() for each in negative_lower]
all_n = negative_upper + negative_cap + negative_lower
punctuation = string.punctuation
punctuation = punctuation.replace('-', '')
print("Podaj tryb działania: ")
print("0. Uruchom opcje 1-4.")
print("1. Stwórz plik, zawierający wszystkie informacje.")
print("2. Stwórz plik zawierający tylko macierz możliwych wariacji polimerów.")
print("3. Podaj liczbę możliwych ułożeń monomerów w polimerze z powtórzeniami (wariacje).")
print("4. Podaj liczbę możliwych unikalnych składów dla polimeru (kombinacje z powtórzeniami).")
print("5. Znajdź polimery mające tą samą kompozycję.")
mode = int(input("Tryb nr: "))
if mode != 5:
    variance = int(input("Podaj ile unikatowych monomerów użyjesz: "))
    size = int(input("Teraz podaj (liczbę) jak duży będzie polimer: "))
else:
    while True:
        podany_polimer = input(
            "Podaj sekwencję polimeru z myślnikami, między monomerami np.: 1-2-3, na którego podstawie chcesz znaleźć "
            "polimery z tym samym składem: ")
        nie_cyfry = punctuation + string.whitespace + string.ascii_letters
        if any(monomer in nie_cyfry for monomer in podany_polimer):
            podany_polimer = podany_polimer.translate({ord(i): None for i in nie_cyfry})
            print("Usuwam interpunkcję, spację oraz litery...")
            response = input(f"Czy sekwencja '{podany_polimer}' jest prawidłowa? ")
            if response in all_p:
                break

        else:
            response = input(f"Czy sekwencja '{podany_polimer}' jest prawidłowa? ")
            if response in all_p:
                break
    size = len(podany_polimer.split('-'))
    podobne = []
    for xs in itertools.product(podany_polimer.split('-'), repeat=size):
        podobne.append(xs)
    print(podobne)
    input('Naciśnij "enter", by zamknąć program.')
    sys.exit("Program został zakończony")
    # response = input(f"Czy do złożenia polimeru mamy '{len(set(podany_polimer.split('-')))}' typów monomerów? ")
    # while True:
    #     if response in all_p:
    #         variance = len(set(podany_polimer.split('-')))
    #         break
    #     else:
    #         variance = input("Podaj ile unikatowych monomerów użyjesz: ")
    #         variance = variance.translate({ord(i): None for i in nie_cyfry})
    #         print(f'Użyto {variance} typów monomerów.')
    #         variance = int(variance)
    #         break
polimer = Inicjalizacja(mode, variance, size)
if mode == 0:
    polimer.variations_with_repeats()
    polimer.combinations_with_repeats()
    nmerlist = NmerList(mode, variance, size).create_nmer_list()
    macierz = Matrix(nmerlist).create_visual_matrix()
    file = File(variance, size, mode, macierz)
    file.create_file_vis()
    file.map_loop_vis()
    file.move_file_vis()
    print('Matryca składu i permutacji została zapisana w: ')
    print('E:\\Bioinformatyka\\PORT\\N-mery\\visual_info\\' + file.file_name)
    macierz2 = Matrix(nmerlist).create_order_matrix()
    file = File(variance, size, mode, macierz2)
    file.create_file_ord()
    file.map_loop_ord()
    file.move_file_ord()
    print('Matryca permutacji została zapisana w: ')
    print('E:\\Bioinformatyka\\PORT\\N-mery\\order_matrix\\' + file.file_name)
    polimer.variations_with_repeats()
    polimer.combinations_with_repeats()
    print(polimer)
elif mode == 1:
    polimer.variations_with_repeats()
    polimer.combinations_with_repeats()
    nmerlist = NmerList(mode, variance, size).create_nmer_list()
    macierz = Matrix(nmerlist).create_visual_matrix()
    file = File(variance, size, mode, macierz)
    file.create_file_vis()
    file.map_loop_vis()
    file.move_file_vis()
    print('Matryca została zapisana w: ')
    print('E:\\Bioinformatyka\\PORT\\N-mery\\visual_info\\' + file.file_name)
elif mode == 2:
    nmerlist = NmerList(mode, variance, size).create_nmer_list()
    macierz = Matrix(nmerlist).create_order_matrix()
    file = File(variance, size, mode, macierz)
    file.create_file_ord()
    file.map_loop_ord()
    file.move_file_ord()
    print('Matryca permutacji została zapisana w: ')
    print('E:\\Bioinformatyka\\PORT\\N-mery\\order_matrix\\' + file.file_name)
elif mode == 3:
    polimer.variations_with_repeats()
    print(polimer)
elif mode == 4:
    polimer.combinations_with_repeats()
    print(polimer)

# kolejny krok
# dodajmy deskryptory polimerów do konwertera

# polarność, kwasowość, wielkość monomerów
# te które najmniej wpływają na widmo
# fluorescencja zależy od środowiska
# dynamika, długie cząsteczki w roztworach
# znana cząsteczka pentamer, sprawdzamy
# mikro zmiany i ogromne zmiany
# zmiany stosunków
