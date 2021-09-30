#!/usr/bin/env python3
# -*- coding: utf-8 -*-

class Point:
    def __init__(self, x, y):
        self.__x = x
        self.__y = y

    def deplace(self, dx, dy):
        self.__x = self.__x + dx
        self.__y = self.__y + dy

    def get_x(self):
        return self.__x

    def get_y(self):
        return self.__y


class PointA(Point):
    def __init__(self, x, y):
        super().__init__(x, y)

    def affiche(self):
        return f"({self.get_x()}, {self.get_y()})"


A = PointA(1, 2)
print(A.affiche())


class Point2:
    def __init__(self, x, y):
        self.__x = x
        self.__y = y

    def aff_coord(self):
        print("Coordonnees : ", self.__x, self.__y)

class PointNom(Point2):
	def __init__(self, nom, x, y):
		super().__init__(x, y)
		self.__nom = nom

	def aff_coord_nom(self):
		print("Nom : ", self.__nom)
		self.aff_coord()

B = PointNom('B', 2, 3)
B.aff_coord_nom()


class Point3:
    def __init__(self, x, y):
        self.__x = x
        self.__y = y

    def affiche(self):
        print("Coordonnees : ", self.__x, self.__y)

class PointNom2(Point3):
	def __init__(self, nom, x, y):
		super().__init__(x, y)
		self.__nom = nom

	def affiche(self):
		print("Nom : ", self.__nom)
		super().affiche()

liste = [PointNom2('A', 0, 0), PointNom2('B', 1, 1), PointNom2('C', 2, 2)]

for point in liste:
	point.affiche()