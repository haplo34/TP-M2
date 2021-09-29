#!/usr/bin/env python3
# -*- coding: utf-8 -*-

class Point:

    def __init__(self, nom, x):
        self.__nom = nom
        self.__x = x
        self.__origine = 0

    def set_origine(self, origine):
        self.__origine = origine

    def get_origine(self):
        return self.__origine
    
    def affiche(self):
        print("Point "
        +self.__nom+" - abscisse = "
        +str(self.__x - self.get_origine())
        +"\nrelative à une origine d'abscisse absolue "
        +str(self.get_origine()))


a = Point('A',3)
a.affiche()
b = Point('B',6)
b.affiche()
a.set_origine(2)
print("On a placé l'origine en", a.get_origine())
a.affiche()
b.affiche()
