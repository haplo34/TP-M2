#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re

class Livre:
    def __init__(self, titre, auteur, nb_pages):
        self.__titre = titre
        self.__auteur = auteur
        self.__nb_pages = nb_pages

    def __str__(self):
        return f"Auteur: {self.__auteur}, " \
            f"Titre: {self.__titre}, " \
            f"{self.__nb_pages} pages."

livre = Livre("Gagner la guerre", "Jean-Philippe Jaworski", 979)
print(livre.__str__())