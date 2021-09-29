#!/usr/bin/env python3
# -*- coding: utf-8 -*-

class Livre:

    def __init__(self, titre, auteur, nb_pages):
        self.__titre = titre
        self.__auteur = auteur
        self.__nb_pages = nb_pages

    def __str__(self):
        return f"Auteur: {self.__auteur}, " \
            f"Titre: {self.__titre}, " \
            f"{self.__nb_pages} pages."

    def get_titre(self):
        return self.__titre

    def get_nb_pages(self):
        return self.__nb_pages

    def set_titre(self, titre):
        self.__titre = titre

    def set_auteur(self, auteur):
        self.__auteur = auteur

    def set_nb_pages(self, nb_pages):
        if nb_pages > 0:
            self.__nb_pages = nb_pages
        else:
            print("Le nombre de pages doit être positif.")

"""
livre_1 = Livre("Gagner la Guerre", "Jean-Philippe Jaworski", 979)
livre_2 = Livre("Dragon d'un Crépuscule d'Automne", "Weis & Hickman", 375)
print(livre_1.__str__())
print(livre_2.__str__())
livre_1.set_nb_pages(10)
livre_2.set_nb_pages(20)
print(livre_1.get_nb_pages())
print(livre_2.get_nb_pages())
print(str(livre_1.get_nb_pages()+livre_2.get_nb_pages()))
"""


class PointAxe:
    def __init__(self, nom, x):
        self.__nom = nom
        self.__x = x

    def affiche(self):
        print("("+self.__nom+", "+str(self.__x)+")")

    def translate(self, dx):
        self.__x += dx        


A = PointAxe('A', 1)
A.affiche()
A.translate(-12)
A.affiche()