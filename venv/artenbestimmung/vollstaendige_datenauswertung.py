import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from scipy.cluster.hierarchy import dendrogram, linkage, complete
import csv
import openpyxl
import ellenberg
import os

class Art:
    def __init__(self, name_l, dataframe=None, name_d=None):
        self.name_l = name_l
        self.name_d = name_d

        self.df = dataframe
        self.stetigkeit = 0

    def bestimmungStetigkeit(self):
        for element in self.df.loc[self.name_l, :].values.tolist():
            if element > 0:
                self.stetigkeit +=1


    def set_df(self, df):
        self.df = df

class Messung:
    def __init__(self, name, dataframe=None):
        self.zeigerwerte = None
        self.name = name

        self.artenzahl = 0
        self.totale_abundanz = 0
        self.shannon_wiener = 0
        self.eveness = 0

        self.array_arten = []

        self.df = dataframe
        self.dominanz_df = None

    def bestimmung_artenzahl_gesamteAbundanz(self):

        for element in self.df[self.name].values.tolist():
            # element entspricht einem prozentualen Wert, weshalb nur Arten die Artenzahl erhöhen, wo auch welche
            # vorhanden waren --> nicht vergessen, self.df enthält alle Messungen
            print(element, type(element))
            if element > 0:
                self.artenzahl +=1
                self.totale_abundanz += element

        # self.array_arten enthält alle Artennamen, die gefunden wurden
        self.array_arten = self.df[self.df[self.name] > 0].index.tolist()

    # bekommt als parameter das Array mit allen Artennamen mitgegeben, weil ich nicht wusste, wie man
    # aus einem dataframe einfach alle zeilen-titel bekommt
    def generiere_dominanz_df(self, arten_namen_array):
        array = self.df[self.name].values.tolist()

        for i in range(0, len(array)):
            # berechnung der dominanz und direktem setzen im array
            array[i] = array[i]/self.totale_abundanz

            # hoffe an dieser stelle, dass zur berechnung des shannon wieners de dominanz verwendet wird
            if array[i] !=0:
                self.shannon_wiener -=array[i] * math.log(array[i])

        # berechnung der eveness
        if self.artenzahl != 0 and self.artenzahl != 1:
            self.eveness = self.shannon_wiener/math.log(self.artenzahl)
        else:
            print("Achtung, Artenzahl = 0 oder 1!")

        # sortieren des dominanz-dataframes
        self.dominanz_df = pd.DataFrame(data=array, columns=[self.name], index=arten_namen_array)
        self.dominanz_df = self.dominanz_df.sort_values(by=[self.name], ascending=False, axis="index")

    def set_df(self, df):
        self.df = df


class allData:
    def __init__(self, pfad_lade_csv, pfad_speichere_excel):
        messungen_namen = []
        ergebnisse = []
        # Transformieren der Braun-Blanquet Werte zu prozentualen Werten der Verteilung
        self.liste_symbole = ["/", "r", "+", "1", "2m", "2a", "2b", "3", "4", "5"]  # Symbole
        self.liste_werte = [0, 0.0005, 0.008, 0.02, 0.04, 0.1, 0.205, 0.38, 0.63, 0.88]  # Prozentuale Werte


        counter = 0
        with open(pfad_lade_csv) as csvfile:
            readCSV = csv.reader(csvfile, delimiter=';')
            next(readCSV, None)
            next(readCSV, None)

            # Start mit Row 3, kein Plan warum
            for row in readCSV:
                # Generierung neues Array für jede Reihe => entspricht ArtenErgebnis
                art_ergebnis_array = []
                # Erste Reihe ist numerierte Messung
                if counter == 0:
                    for i in range(3, len(row)):
                        if row[i] != "":
                            messungen_namen.append(int(row[i]))
                    counter +=1
                else:
                    # Nun nehmen wir die Arten, der name_l am Index 0
                    art_ergebnis_array.append(row[0])
                    # die lablanque-Werte stehen ab Zeile 3
                    for i in range(3, len(row)):
                        if row[i] != "":
                            art_ergebnis_array.append(row[i])
                    ergebnisse.append(art_ergebnis_array)

        # Weil Arrays toll sind
        self.liste_messungen = []
        self.liste_arten = []
        self.liste_artenname = []

        # Aus Ergebnissen werden neue Arten generiert und liste_artenname/liste_arten hinzugefügt
        # Artenname wird aus Ergebnis Array entfernt, aber warum?
        for element in ergebnisse:
            neueArt = Art(element[0])
            self.liste_artenname.append(element[0])
            self.liste_arten.append(neueArt)
            element.pop(0)

        # Generierung Messungs-Objekt
        for element in messungen_namen:
            neueMessung = Messung(name=element)
            self.liste_messungen.append(neueMessung)

        # Dataframe mit Symbolen, falls man die Ausdrucken will
        self.data_frame_symbole = pd.DataFrame(data=ergebnisse, columns=messungen_namen, index=self.liste_artenname)
        # Oben wurde der Datenname entfernt, damit man hier nun im Array Ergebnisse einfacher Blanque zu Zahlen
        for element in ergebnisse:
            for i in range(0, len(element)):
                for counter, value in enumerate(self.liste_symbole):
                    if element[i] == value:
                        element[i] = self.liste_werte[counter]

        # Plott in proz. Zahlenwerten
        self.df = pd.DataFrame(data=ergebnisse, columns=messungen_namen, index=self.liste_artenname)
        print(self.df)
        # Jedes "Messungs"-Objekt bekommt Property des proz. Dataframes
        for element in self.liste_messungen:
            element.set_df(self.df)

            # Jede Messung bestimmt ihre Artenzahl und Abundanz
            element.bestimmung_artenzahl_gesamteAbundanz()
            # Und generiert ihre Dominanz_dataframe und verwendet hierzu als parameter die liste aller Artennamen
            element.generiere_dominanz_df(self.liste_artenname)

        liste_stetigkeit = []
        for element in self.liste_arten:
            element.set_df(self.df)
            element.bestimmungStetigkeit()
            liste_stetigkeit.append(element.stetigkeit)
        self.df["Stetigkeit"] = liste_stetigkeit
        self.df = self.df.sort_values(by=["Stetigkeit"], ascending=False, axis="index")

        self.jaccard_df = self.generierung_jaccard_df()

        self.ruzicka_df = self.generierung_ruzicka()
        self.generate_zeigerwerte()

        self.generate_excel(pfad_speichere_excel)


    def generate_zeigerwerte(self):
        for element in self.liste_messungen:
            path = "zeigerwerte/Messung_"+str(element.name)

            element.zeigerwerte = ellenberg.do_all(element.array_arten, path)

    def generierung_jaccard_df(self):
        new_df = pd.DataFrame(index=[element.name for element in self.liste_messungen],
                              columns=[element.name for element in self.liste_messungen])
        for row in self.liste_messungen:
            for column in self.liste_messungen:
                if row.name == column.name:
                    new_df._set_value(row.name, column.name, "/")
                else:
                    anzahl_gemeinsame_elemente = len(set(row.array_arten).intersection(column.array_arten))
                    anzahl_eigene_elemente_row = len([element for element in row.array_arten if element not in column.array_arten])
                    anzahl_eigene_elemente_column = len([element for element in column.array_arten if element not in row.array_arten])

                    jaccard = round(100*anzahl_gemeinsame_elemente/(anzahl_eigene_elemente_column+anzahl_gemeinsame_elemente+anzahl_eigene_elemente_row), 2)
                    new_df._set_value(row.name, column.name, jaccard)

        return new_df

    def generierung_ruzicka(self):
        new_df = pd.DataFrame(index=[element.name for element in self.liste_messungen],
                              columns=[element.name for element in self.liste_messungen])

        for row in self.liste_messungen:

            for column in self.liste_messungen:

                summe_zähler = 0
                summe_nenner = 0

                if row.name == column.name:
                    new_df._set_value(row.name, column.name, "/")

                else:
                    row_array = self.df[row.name].values.tolist()
                    column_array = self.df[column.name].values.tolist()

                    for counter, value in enumerate(row_array):
                        if value >= column_array[counter]:
                            summe_nenner += value
                            summe_zähler += column_array[counter]
                        else:
                            summe_nenner += column_array[counter]
                            summe_zähler += value

                        if 100*summe_zähler/summe_nenner == None:
                            print(summe_zähler, summe_nenner)
                    new_df._set_value(row.name, column.name, round(100*summe_zähler/summe_nenner, 2))

        return new_df


    def generate_excel(self, path):

        writer = pd.ExcelWriter(path+".xlsx")
        book = writer.book

        self.df.to_excel(writer, sheet_name="Alle Messungen")
        self.jaccard_df.to_excel(writer, sheet_name="Jaccard-Matrix")
        self.ruzicka_df.to_excel(writer, sheet_name="Ruzicka-Matrix")

        array_zeigerwerte = ['Licht: ', 'Temperatur: ', 'Kontinentalitaet: ', 'Feuchte: ', 'Reaktion: ', 'Stickstoff: ']

        for element in self.liste_messungen:

            data_array = ["Artenzahl: "+str(element.artenzahl), "Shannon-Wiener: "+str(element.shannon_wiener), "Evenness: "+str(element.eveness)]

            for counter, value in enumerate(array_zeigerwerte):
                data_array.append("Mean Zeigerwert " + value + str(element.zeigerwerte[counter]["mean"]) + "-" + "Standardabweichung: "+ str(element.zeigerwerte[counter]["dev"])+ "Anzahl n: "+ str(element.zeigerwerte[counter]["count"]))


            messung = pd.DataFrame(data=data_array, columns=["Allgemeine Daten:"])

            sheet_name = "Messung " + str(element.name) + " Dominanz-Rang"
            messung.to_excel(writer, sheet_name=sheet_name)
            element.dominanz_df.to_excel(writer, sheet_name=sheet_name, startrow=12)

            ws = book.get_sheet_by_name(sheet_name)


            copy_array = element.dominanz_df[element.name].values.tolist()

            kumulative_werte = []
            index = []
            kumulativ = 0
            counter = 1
            for i in range(0, len(copy_array)):
                if copy_array[i] != 0:
                    kumulativ += copy_array[i]
                    kumulative_werte.append(kumulativ*100)
                    index.append(counter)
                    counter +=1

            plt.figure()
            plt.plot(index, kumulative_werte, "-ro")
            plt.xticks(index)
            plt.ylabel("Kumulative Dominanz [%]")
            plt.xlabel("Rang")
            plot_name = "Plot k-Dominanz-Rang - Messung " + str(element.name)
            plt.title(plot_name)
            plt.savefig(plot_name, dpi=100)
            img = openpyxl.drawing.image.Image(plot_name+".png")

            ws.add_image(img, "D6")

        ws = book.get_sheet_by_name("Jaccard-Matrix")
        path_jaccard = "Jaccard"
        self.create_dendrogram(self.jaccard_df, path_jaccard)
        img = openpyxl.drawing.image.Image(path_jaccard+".png")
        ws.add_image(img, "A9")

        ws = book.get_sheet_by_name("Ruzicka-Matrix")
        path_ruzicka = "Ruzicka"
        self.create_dendrogram(self.ruzicka_df, path_ruzicka)
        img = openpyxl.drawing.image.Image(path_ruzicka+".png")
        ws.add_image(img, "A9")

        writer.save()

    def create_dendrogram(self, df, plot_name):
        matrix = df.to_numpy()
        matrix_without = []

        for counter, element in enumerate(matrix):
            for count, val in enumerate(element):
                if val == "/":
                    pass
                elif counter < count:
                    matrix_without.append(val)

        for i in range(0, len(matrix_without)):
            matrix_without[i] = 100 - matrix_without[i]


        Z = linkage(matrix_without, "complete", optimal_ordering=True)
        # Z = complete(matrix_without)

        liste_labels = ["Prob. " + str(messung.name) for messung in self.liste_messungen]
        #liste_labels = ["eG1", "eG2", "eG3", "eG4", "B5", "B6", "W7", "W8", "W9", "W10", "nW11", "nW12"]
        fig = plt.figure()

        dn = dendrogram(Z, orientation="right", labels=liste_labels)

        plt.xlim(0, 105)
        #plt.xlabel = (plot_name+" Ähnlichkeit [%]")
        plt.gca().set_xticklabels([100, 80, 60, 40, 20, 0])
        plt.gca().set_xlabel(plot_name+" Ähnlichkeit [%]")

        plt.savefig(plot_name, dpi=150)

#dirname = os.path.dirname(__file__)
#dirIn = dirname +"/Daten/Input_fuerDatenauswertung.csv"
#dirOut = dirname+"Daten/generierteAuswertung"
lets_do_this = allData(pfad_lade_csv="/Users/leonkollarczik/Desktop/artenbestimmung/venv/Daten/Input_fuerDatenauswertung.csv",
               pfad_speichere_excel="/Users/leonkollarczik/Desktop/artenbestimmung/venv/Daten/generierteAuswertung")
