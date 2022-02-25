# -*- coding: utf-8 -*-
# units.py by Gustav Lindberg
# https://github.com/Gustav-Lindberg/Units.py

from fractions import Fraction as _Fraction
from unicodedata import lookup as _unicode
from numpy.core.numeric import isclose as _isclose
from numpy import pi
from numpy import exp as _exp, log as _log, log10 as _log10
from numpy import sin as _sin, cos as _cos, tan as _tan
from numpy import arcsin as _arcsin, arccos as _arccos, arctan as _arctan
import re as _re

class Unit:
    def __init__(this, other=None):
        try:
            if(isinstance(other, Quantity)):
                this._kilograms = other._unit._kilograms
                this._meters = other._unit._meters
                this._seconds = other._unit._seconds
                this._amperes = other._unit._amperes
                this._kelvins = other._unit._kelvins
                this._moles = other._unit._moles
                this._multiple = other._number * other._unit._multiple
                return
        except NameError:
            pass
        
        if(isinstance(other, int) or isinstance(other, float) or isinstance(other, complex)):
            multiple = other
            other = None
        elif(other == None):
            multiple = 1.0
        
        if(other != None and not(isinstance(other, Unit))):
            raise TypeError("Cannot convert {} to unit".format(type(other)))
            
        this._kilograms = 0 if other == None else other._kilograms
        this._meters = 0 if other == None else other._meters
        this._seconds = 0 if other == None else other._seconds
        this._amperes = 0 if other == None else other._amperes
        this._kelvins = 0 if other == None else other._kelvins
        this._moles = 0 if other == None else other._moles
        this._multiple = multiple if other == None else other._multiple
        
    def __copy__(this):
        return Unit(this)
    
    def __deepcopy__(this):
        return Unit(this)
        
        
    def __imul__(this, other):
        if(isinstance(other, Unit) or isinstance(other, Quantity)):
            o = Unit(other)
            this._kilograms += o._kilograms
            this._meters += o._meters
            this._seconds += o._seconds
            this._amperes += o._amperes
            this._kelvins += o._kelvins
            this._moles += o._moles
            this._multiple *= o._multiple
        else:
            this._multiple *= float(other)
        return this
            
            
    def __mul__(this, other):
        try:
            iter(other)
            return [o * this for o in other]
        except Exception:
            pass
        toReturn = Unit(this)
        toReturn *= other
        if(isinstance(other, Unit)):
            return toReturn
        elif(isinstance(other, Quantity)):
            return Quantity(1, toReturn).toUnit(this * other._unit)
        else:
            return Quantity(1, toReturn).toUnit(this)
    
    def __rmul__(this, other):
        return this * other
    
    
    def __ipow__(this, other):
        o = _Fraction(other)
        this._kilograms *= o
        this._meters *= o
        this._seconds *= o
        this._amperes *= o
        this._kelvins *= o
        this._moles *= o
    
        this._multiple **= o
        return this
        
        
    def __pow__(this, other):
        toReturn = Unit(this)
        toReturn **= other
        return toReturn
    
    
    def __idiv__(this, other):
        this *= other ** -1
        return this
        
    def __truediv__(this, other):
        try:
            iter(other)
            return [o / this for o in other]
        except Exception:
            pass
        return this * other ** -1
    
    def __rtruediv__(this, other):
        return other * this ** -1
    
    
    def isMultiple(this, other):
        return (
            isinstance(other, Unit) and
            this._kilograms == other._kilograms and
            this._meters == other._meters and
            this._seconds == other._seconds and
            this._amperes == other._amperes and
            this._kelvins == other._kelvins and
            this._moles == other._moles
        )
        
    def isSiUnit(this):
        return _isclose(this._multiple, 1.0)
    
    def isPlanckUnit(this):
        return Quantity(1, this).toPlanckUnits()._unit == this
            
    def __eq__(this, other):
        if(isinstance(other, Quantity)):
            return this == Unit(other)
        return this.isMultiple(other) and _isclose(this._multiple / other._multiple, 1)
        
    def __ne__(this, other):
        return not(this == other)
    
    def __nonzero__(this):
        return this._multiple != 0
    
    
    def __neg__(this):
        return this * -1
    
    def __pos__(this):
        return this
    
    
    def __hash__(this):
        return hash((
            this._kilograms,
            this._meters,
            this._seconds,
            this._amperes,
            this._kelvins,
            this._moles,
            this._multiple
        ))
        
    def __int__(this):
        return int(float(this))
        
    def __float__(this):
        if(not(this.isMultiple(dimensionless))):
            raise ValueError("Only units that are multiples of dimensionless can be converted to float or complex, not " + str(this))
        return this._multiple
    
    def __repr__(this):
        return str(this)
    
    def __str__(this):
        if(this == dimensionless):
            return "1"
        elif(this._moles == 1 and this != mol):
            thisTimesNA = str(Unit(this * NA))
            if(not thisTimesNA[0].isdigit() or thisTimesNA == "1"):
                return thisTimesNA
            return str(this / mol) + "·mol"
        elif(this._moles == -1 and this != mol):
            thisDividedByNA = str(Unit(this / NA))
            if(not thisDividedByNA[0].isdigit() or thisDividedByNA == "1"):
                return thisDividedByNA
            toReturn = str(this * mol)
            if('/' in toReturn):
                toReturn += "/mol"
            else:
                toReturn += "·mol" + _unicode("SUPERSCRIPT MINUS") + _unicode("SUPERSCRIPT ONE")
            return toReturn
        elif(this == t):
            return "t"
        elif(this == minute):
            return "min"
        elif(this == minute**-1):
            return "min" + _unicode("SUPERSCRIPT MINUS") + _unicode("SUPERSCRIPT ONE")
        elif(this == hour):
            return "hr"
        elif(this == hour**-1):
            return "hr" + _unicode("SUPERSCRIPT MINUS") + _unicode("SUPERSCRIPT ONE")
        elif(this == d):
            return "d"
        elif(this == d**-1):
            return "d" + _unicode("SUPERSCRIPT MINUS") + _unicode("SUPERSCRIPT ONE")
        elif(this == w):
            return "w"
        elif(this == w**-1):
            return "w" + _unicode("SUPERSCRIPT MINUS") + _unicode("SUPERSCRIPT ONE")
        elif(this == y):
            return "y"
        elif(this == y**-1):
            return "y" + _unicode("SUPERSCRIPT MINUS") + _unicode("SUPERSCRIPT ONE")
        elif(this == Å):
            return "Å"
        elif(this == Å**2):
            return "Å²"
        elif(this == Å**3):
            return "Å³"
        elif(this == au):
            return "au"
        elif(this == au**2):
            return "au²"
        elif(this == au**3):
            return "au³"
        elif(this == ly):
            return "ly"
        elif(this == ly**2):
            return "ly²"
        elif(this == ly**3):
            return "ly³"
        elif(this == ha):
            return "ha"
        elif(this == b):
            return "b"
        elif(this == inch):
            return "in"
        elif(this == inch**2):
            return "in²"
        elif(this == inch**3):
            return "in³"
        elif(this == ft):
            return "ft"
        elif(this == ft**2):
            return "ft²"
        elif(this == ft**3):
            return "ft³"
        elif(this == yd):
            return "yd"
        elif(this == yd**2):
            return "yd²"
        elif(this == yd**3):
            return "yd³"
        elif(this == mi):
            return "mi"
        elif(this == mi**2):
            return "mi²"
        elif(this == mi**3):
            return "mi³"
        elif(this == tsp):
            return "tsp"
        elif(this == tbsp):
            return "tbsp"
        elif(this == floz):
            return "fl oz"
        elif(this == cp):
            return "cp"
        elif(this == pt):
            return "pt"
        elif(this == qt):
            return "qt"
        elif(this == oz):
            return "oz"
        elif(this == lb):
            return "lb"
        elif(this == u):
            return "u"
        elif(this == km / hour):
            return "km/h"
        elif(this == mi / hour):
            return "mph"
        elif(this == kW * hour):
            return "kWh"
        elif(this == atm):
            return "atm"
        elif(this == degree):
            return _unicode("DEGREE SIGN")
        elif(this == Unit(e)):
            return "e"
        elif(this == Unit(G)):
            return "G"
        elif(this == Unit(hbar)):
            return _unicode("LATIN SMALL LETTER H WITH STROKE")
        elif(this == Unit(c)):
            return "c"
        elif(this == Unit(epsilon0)):
            return _unicode("GREEK SMALL LETTER EPSILON") + _unicode("SUBSCRIPT ZERO")
        elif(this == Unit(mu0)):
            return _unicode("GREEK SMALL LETTER MU") + _unicode("SUBSCRIPT ZERO")
        elif(this == Unit(planckEnergy)):
            return "E" + str(_unicode("LATIN SUBSCRIPT SMALL LETTER P"))
        unitToString = dict()
        unitToString[m] = "m"
        unitToString[s] = "s"
        unitToString[A] = "A"
        unitToString[K] = "K"
        unitToString[Hz] = "Hz"
        unitToString[N] = "N"
        unitToString[Pa] = "Pa"
        unitToString[J] = "J"
        unitToString[W] = "W"
        unitToString[C] = "C"
        unitToString[V] = "V"
        unitToString[F] = "F"
        unitToString[Ohm] = _unicode("GREEK CAPITAL LETTER OMEGA")
        unitToString[S] = "S"
        unitToString[Wb] = "Wb"
        unitToString[T] = "T"
        unitToString[H] = "H"
        basicUnits = dict(unitToString)
        unitToString[L] = "L"
        unitToString[g] = "g"
        unitToString[eV] = "eV"
        unitToString[Unit(eV / c)] = "eV/c"
        unitToString[Unit(eV / c**2)] = "eV/c" + _unicode("SUPERSCRIPT TWO")
        unitToString[mol] = "mol"
        unitToString[Da] = "Da"
        unitToString[cal] = "cal"
        unitToString[pc] = "pc"
        for unit in list(basicUnits.keys()) + [g, L, Unit(eV), Unit(eV / c), Unit(eV / c**2), mol, Da, cal, pc]:
            unitToString[Unit(1e-1 * unit)] = "d" + unitToString[unit]
            unitToString[Unit(1e-2 * unit)] = "c" + unitToString[unit]
            unitToString[Unit(1e-3 * unit)] = "m" + unitToString[unit]
            unitToString[Unit(1e-6 * unit)] = _unicode("GREEK SMALL LETTER MU") + unitToString[unit]
            unitToString[Unit(1e-9 * unit)] = "n" + unitToString[unit]
            unitToString[Unit(1e-12 * unit)] = "p" + unitToString[unit]
            unitToString[Unit(1e-15 * unit)] = "f" + unitToString[unit]
            unitToString[Unit(1e-18 * unit)] = "a" + unitToString[unit]
            unitToString[Unit(1e-21 * unit)] = "z" + unitToString[unit]
            unitToString[Unit(1e-24 * unit)] = "y" + unitToString[unit]
            unitToString[Unit(1e1 * unit)] = "da" + unitToString[unit]
            unitToString[Unit(1e2 * unit)] = "h" + unitToString[unit]
            if(unit != L):
                unitToString[Unit(1e3 * unit)] = "k" + unitToString[unit]
            unitToString[Unit(1e6 * unit)] = "M" + unitToString[unit]
            unitToString[Unit(1e9 * unit)] = "G" + unitToString[unit]
            unitToString[Unit(1e12 * unit)] = "T" + unitToString[unit]
            unitToString[Unit(1e15 * unit)] = "P" + unitToString[unit]
            unitToString[Unit(1e18 * unit)] = "E" + unitToString[unit]
            unitToString[Unit(1e21 * unit)] = "Z" + unitToString[unit]
            unitToString[Unit(1e24 * unit)] = "Y" + unitToString[unit]
        basicUnits[kg] = "kg"
        
        def toSuperscript(string):
            toReturn = str(string)
            names = ["ZERO", "ONE", "TWO", "THREE", "FOUR", "FIVE", "SIX", "SEVEN", "EIGHT", "NINE"]
            for i in range(10):
                toReturn = toReturn.replace(str(i), _unicode("SUPERSCRIPT {}".format(names[i])))
            toReturn = toReturn.replace('-', _unicode("SUPERSCRIPT MINUS"))
            toReturn = toReturn.replace('/', _unicode("RIGHT RAISED OMISSION BRACKET"))
            return toReturn
        
        try:
            if(this.isPlanckUnit()):
                toReturn = ""
                if(this._kilograms == 1):
                    toReturn += "m{}·".format(_unicode("LATIN SUBSCRIPT SMALL LETTER P"))
                elif(this._kilograms != 0):
                    toReturn += "m{}{}·".format(_unicode("LATIN SUBSCRIPT SMALL LETTER P"), toSuperscript(this._kilograms))
                if(this._meters == 1):
                    toReturn += "l{}·".format(_unicode("LATIN SUBSCRIPT SMALL LETTER P"))
                elif(this._meters != 0):
                    toReturn += "l{}{}·".format(_unicode("LATIN SUBSCRIPT SMALL LETTER P"), toSuperscript(this._meters))
                if(this._seconds - this._amperes == 1):
                    toReturn += "t{}·".format(_unicode("LATIN SUBSCRIPT SMALL LETTER P"))
                elif(this._seconds - this._amperes != 0):
                    toReturn += "t{}{}·".format(_unicode("LATIN SUBSCRIPT SMALL LETTER P"), toSuperscript(this._seconds - this._amperes))
                if(this._amperes == 1):
                    toReturn += "q{}·".format(_unicode("LATIN SUBSCRIPT SMALL LETTER P"))
                elif(this._amperes != 0):
                    toReturn += "q{}{}·".format(_unicode("LATIN SUBSCRIPT SMALL LETTER P"), toSuperscript(this._amperes))
                if(this._kelvins == 1):
                    toReturn += "T{}·".format(_unicode("LATIN SUBSCRIPT SMALL LETTER P"))
                elif(this._kelvins != 0):
                    toReturn += "T{}{}·".format(_unicode("LATIN SUBSCRIPT SMALL LETTER P"), toSuperscript(this._kelvins))
                if(this._moles != 0 and toReturn == ""):
                    return "1"
                return toReturn[:-1]
        except OverflowError:
            pass
            
        if(this in unitToString.keys()):
            return unitToString[this]
        for unit in basicUnits.keys():
            if(this.isMultiple(unit)):
                return str(this._multiple) + unitToString[unit]
        if(this == dimensionless):
            return "1"
        if(this.isMultiple(dimensionless)):
            return str(this._multiple)
        for unit in unitToString.keys():
            for i in range(-3, 4):
                if(this == unit ** i):
                    return unitToString[unit] + toSuperscript(i)
        for unit in basicUnits.keys():
            for i in range(-3, 4):
                if(this.isMultiple(unit ** i)):
                    return str(this._multiple) + unitToString[unit] + toSuperscript(i)
        units = [m, s, J, eV, V, K, F, H, W, C]
        for unit1 in units:
            for unit2 in units:
                for i in [1, 2, 3]:
                    if(this == unit1 * unit2**i):
                        return unitToString[unit1] + "·" + unitToString[unit2] + ("" if i == 1 else toSuperscript(i))
                    elif(this == unit1 / unit2**i):
                        return unitToString[unit1] + "/" + unitToString[unit2] + ("" if i == 1 else toSuperscript(i))
        
        toReturn = "" if this.isSiUnit() else str(this._multiple)
        if(this._kilograms == 1):
            toReturn += "kg·"
        elif(this._kilograms != 0):
            toReturn += "kg" + toSuperscript(this._kilograms) + "·"
        if(this._meters == 1):
            toReturn += "m·"
        elif(this._meters != 0):
            toReturn += "m" + toSuperscript(this._meters) + "·"
        if(this._seconds == 1):
            toReturn += "s·"
        elif(this._seconds != 0):
            toReturn += "s" + toSuperscript(this._seconds) + "·"
        if(this._amperes == 1):
            toReturn += "A·"
        elif(this._amperes != 0):
            toReturn += "A" + toSuperscript(this._amperes) + "·"
        if(this._kelvins == 1):
            toReturn += "K·"
        elif(this._kelvins != 0):
            toReturn += "K" + toSuperscript(this._kelvins) + "·"
        if(this._moles == 1):
            toReturn += "mol·"
        elif(this._moles != 0):
            toReturn += "mol" + toSuperscript(this._moles) + "·"
        return toReturn[:-1]
    
class Quantity:
    def __init__(this, number=1.0, unit=Unit()):
        this._unit = Unit(unit)
        if(isinstance(number, Quantity)):
            this._unit *= number._unit
            this._number = number._number
        elif(isinstance(number, Unit)):
            this._unit = Unit(number)
            this._number = 1
        elif(isinstance(number, complex)):
            this._number = number
        else:
            this._number = float(number)
    
    def __copy__(this):
        toReturn = Quantity
        toReturn._number = this._number
        toReturn._unit = this._unit
        return toReturn
    
    def __deepcopy__(this):
        return Quantity(this)
        
        
    def __imul__(this, other):
        if(isinstance(other, Quantity)):
            this._number *= other._number
            this._unit *= other._unit
        elif(isinstance(other, Unit)):
            this._unit *= other
        elif(isinstance(other, complex)):
            this._number *= other
        else:
            this._number *= float(other)
        return this
    
    def __mul__(this, other):
        try:
            iter(other)
            return [o * this for o in other]
        except Exception:
            pass
        toReturn = Quantity(this)
        toReturn *= other
        return toReturn.dimensionlessToFloat()
    
    def __rmul__(this, other):
        return this * other
    
    def __ipow__(this, other):
        this._number **= other
        this._unit **= other
        return this
        
    def __pow__(this, other):
        try:
            iter(other)
            return [o ** this for o in other]
        except Exception:
            pass
        toReturn = Quantity(this)
        toReturn **= other
        return toReturn.dimensionlessToFloat()
    
    def __rpow__(this, other):
        return other ** float(this)
    
    def __idiv__(this, other):
        this *= other ** -1
        return this
    
    def __truediv__(this, other):
        return this * other ** -1
    
    def __rtruediv__(this, other):
        return other * this ** -1
    
    def __neg__(this):
        return (this * -1).dimensionlessToFloat()
    
    def __pos__(this):
        return this.dimensionlessToFloat()
    
    
    def __iadd__(this, other):
        if(other == 0 and Quantity(other)._unit.isMultiple(dimensionless)):
            return this
        if(this == 0 and Quantity(this)._unit.isMultiple(dimensionless)):
            this._number = other._number
            this._unit = other._unit
            return other
        o = Quantity(other)
        if(not(this._unit.isMultiple(o._unit))):
            raise ValueError("Cannot add two quantities with different units: {} and {}".format(this._unit, o._unit))
        this._number += o._number * float(o._unit / this._unit)
        return this
        
    def __add__(this, other):
        try:
            iter(other)
            return [o + this for o in other]
        except Exception:
            pass
        toReturn = Quantity(this)
        toReturn += other
        return toReturn.dimensionlessToFloat()
    
    def __radd__(this, other):
        return this + other
    
    def __isub__(this, other):
        this += (-other)
        return this
    
    def __sub__(this, other):
        try:
            iter(other)
            return [o - this for o in other]
        except Exception:
            pass
        return this + (-other)
    
    def __rsub__(this, other):
        return -(this - other)
    
    def __eq__(this, other):
        if(other == 0):
            return this._number == 0 or this._unit._multiple == 0
        elif(isinstance(other, Quantity)):
            return this._unit.isMultiple(other._unit) and _isclose((this._number * this._unit._multiple) / (other._number * other._unit._multiple), 1)
        elif(isinstance(other, Unit)):
            return this == Quantity(1, other)
        elif(this._unit.isMultiple(dimensionless) and isinstance(other, complex)):
            return complex(this) == other
        else:
            try:
                return float(this) == float(other)
            except Exception:
                return False
        
    def __neq__(this, other):
        return not(this == other)
    
    def __nonzero__(this):
        return this != 0
    
    def __lt__(this, other):
        if(this == 0 and Quantity(this)._unit.isMultiple(dimensionless)):
            return 0 <= other._number
        if(other == 0 and Quantity(other)._unit.isMultiple(dimensionless)):
            return this._number <= 0
        o = Quantity(other)
        if(not(this._unit.isMultiple(o._unit))):
            raise ValueError("Cannot compare two quantities with different units: {} and {}".format(this._unit, o._unit))
        return this._number * this._unit._multiple < o._number * o._unit._multiple
    
    def __le__(this, other):
        return this < other or this == other
    
    def __gt__(this, other):
        return not(this <= other)
    
    def __ge__(this, other):
        return not(this < other)
    
    def __hash__(this):
        return hash((this._number, this._unit))
    
    def __int__(this):
        return int(complex(this))
        
    def __float__(this):
        return float(this._number) * float(this._unit)
    
    def __complex__(this):
        return complex(this._number) * float(this._unit)
    
    def __repr__(this):
        return str(this)
    
    def __str__(this):
        if(this._unit == dimensionless):
            return str(this._number)
        elif(not(str(this._unit)[0].isdigit()) and str(this._unit)[0] != '.' and str(this._unit)[0] != '-'):
            return str(this._number) + _re.sub("^[0-9.-]+", "", str(this._unit))
        elif(this._unit.isMultiple(dimensionless)):
            return str(this._number * float(this._unit))
        else:
            return str(this._number * this._unit._multiple) + _re.sub("^[0-9.-]+", "", str(Unit(this._unit / this._unit._multiple)))
        
    def toUnit(this, unit):
        unit = Unit(unit)
        if(this._unit._kilograms != unit._kilograms or this._unit._meters != unit._meters or this._unit._seconds != unit._seconds or this._unit._amperes != unit._amperes or this._unit._kelvins != unit._kelvins):
            raise ValueError("Cannot convert quantity in unit {} to quantity in unit {}".format(this._unit, unit))
        toReturn = Quantity(this._number * this._unit._multiple / unit._multiple, unit)
        if(this._unit._moles != unit._moles):
            toReturn *= (NA * mol) ** (this._unit._moles - unit._moles)
        return toReturn
    
    def toSiUnits(this):
        unit = Unit(this._unit)
        unit._multiple = 1
        return this.toUnit(unit)
    
    def toPlanckUnits(this):
        unit = Unit(planckMass ** this._unit._kilograms * planckLength ** this._unit._meters * planckTime ** this._unit._seconds * (planckCharge / planckTime) ** this._unit._amperes * planckTemperature ** this._unit._kelvins)
        return this.toUnit(unit)
    
    def molesToDimensionless(this):
        unit = Unit(this._unit)
        unit._moles = 0
        return this.toUnit(unit)
    
    def isDimensionless(this):
        return this._unit.isMultiple(dimensionless)
    
    def dimensionlessToFloat(this):
        if(this.isDimensionless()):
            try:
                return float(this)
            except TypeError:
                return complex(this)
        else:
            return this
    
    
    def exp(this):
        return _exp(float(this))
    
    def log(this):
        return _log(float(this))
    
    def log10(this):
        return _log10(float(this))
    
    def sin(this):
        return _sin(float(this))
    
    def cos(this):
        return _cos(float(this))
    
    def tan(this):
        return _tan(float(this))
    
    def arcsin(this):
        return _arcsin(float(this))
    
    def arccos(this):
        return _arccos(float(this))
    
    def arctan(this):
        return _arctan(float(this))
    
    def sqrt(this):
        return this ** _Fraction(1, 2)
    
    def __abs__(this):
        toReturn = Quantity(this)
        toReturn._number = abs(toReturn._number)
        return toReturn
    
    def __round__(this):
        return Quantity(round(this._number), this._unit)

dimensionless = Unit()
kg = Unit()
kg._kilograms = 1
m = Unit()
m._meters = 1
s = Unit()
s._seconds = 1
A = Unit()
A._amperes = 1
K = Unit()
K._kelvins = 1
mol = Unit()
mol._moles = 1

Hz = s**-1
N = kg * m * s**-2
Pa = N / m**2
J = N * m
W = J / s
C = A * s
V = W / A
F = C / V
Ohm = V / A
S = Ohm**-1
Wb = V * s
T = Wb / m**2
H = Wb / A
Bq = Hz
Gy = J / kg
Sv = Gy
degree = degrees = Unit(pi / 180)
arcmin = Unit(degree / 60)
arcsec = Unit(degree / 3600)
radian = dimensionless

cm = Unit(1e-2 * m)
mm = Unit(1e-3 * m)
µm = micrometer = Unit(1e-6 * m)
nm = Unit(1e-9 * m)
pm = Unit(1e-12 * m)
fm = Unit(1e-15 * m)
km = Unit(1e3 * m)

MPa = Unit(1e6 * Pa)
kPa = Unit(1e3 * Pa)
MPa = Unit(1e6 * Pa)
GPa = Unit(1e9 * Pa)
TPa = Unit(1e12 * Pa)
mPa = Unit(1e-3 * Pa)
nPa = Unit(1e-9 * Pa)

TJ = Unit(1e12 * J)
GJ = Unit(1e9 * J)
MJ = Unit(1e6 * J)
kJ = Unit(1e3 * J)
mJ = Unit(1e-3 * J)

TW = Unit(1e12 * W)
GW = Unit(1e9 * W)
MW = Unit(1e6 * W)
kW = Unit(1e3 * W)
mW = Unit(1e-3 * W)

GV = Unit(1e9 * V)
MV = Unit(1e6 * V)
kV = Unit(1e3 * V)
mV = Unit(1e-3 * V)

kA = Unit(1e3 * A)
mA = Unit(1e-3 * A)
µA = microAmpere = Unit(1e-6 * A)
nA = Unit(1e-9 * A)
pA = Unit(1e-12 * A)
fA = Unit(1e-15 * A)

mF = Unit(1e-3 * F)
µF = microFarad = Unit(1e-6 * F)
nF = Unit(1e-9 * F)
pF = Unit(1e-12 * F)

mOhm = milliOhm = Unit(1e-3 * Ohm)
kOhm = kiloOhm = Unit(1e3 * Ohm)
MOhm = megaOhm = Unit(1e6 * Ohm)

ms = Unit(1e-3 * s)
µs = microsecond = Unit(1e-6 * s)
ns = Unit(1e-9 * s)
ps = Unit(1e-12 * s)
fs = Unit(1e-15 * s)

mg = Unit(1e-6 * kg)
µg = microgram = Unit(1e-9 * kg)
ng = Unit(1e-12 * kg)

kBq = Unit(1e3 * Bq)
MBq = Unit(1e6 * Bq)
kHz = Unit(1e3 * Hz)
MHz = Unit(1e6 * Hz)
GHz = Unit(1e9 * Hz)
THz = Unit(1e12 * Hz)

mC = Unit(1e-3 * C)
µC = microCoulomb = Unit(1e-6 * C)
nC = Unit(1e-9 * C)
pC = Unit(1e-12 * C)
fC = Unit(1e-15 * C)

kB = boltzmannConstant = Quantity(1.380649e-23, J / K)
e = elementaryCharge = Quantity(1.602176634e-19, C)
G = gravitationalConstant = Quantity(6.67430e-11, m**3 * kg**-1 * s**-2)
h = planckConstant = Quantity(6.62607015e-34, J * s)
hbar = reducedPlanckConstant = h / (2 * pi)
c = speedOfLight = Quantity(299792458.0, m / s)
epsilon0 = Quantity(8.8541878128e-12, F / m)
mu0 = 1 / (epsilon0 * c**2)
NA = avogadroNumber = Quantity(6.02214076e23, mol ** -1)
R = kB * NA
KJ = josephsonConstant = 2 * e / h
RK = vonKlitzingConstant = h / e**2
magneticFluxQuantum = 1 / josephsonConstant
stefanBoltzmannConstant = Quantity(5.670373e-8, W / (m**2 * K**4))

fineStructureConstant = alphaElectromagnetic = 7.2973525693e-3
alphaStrong = 0.1179
alphaWeak = 1/30
gElectromagnetic = (4 * pi * alphaElectromagnetic) ** (1/2)
gStrong = (4 * pi * alphaStrong) ** (1/2)
gWeak = (4 * pi * alphaWeak) ** (1/2)

cabibboAngle = Quantity(13.02, degree)
weakMixingAngle = Quantity(28.2, degree)
electronMass = me = Quantity(9.1093837015e-31, kg)
protonMass = mp = Quantity(1.67262192369e-27, kg)
neutronMass = mn = Quantity(1.674927351e-27, kg)
cosmologicalConstant = Quantity(1.1056e-52, m**-2)
muB = bohrMagneton = e * hbar / (2 * me)
rydbergConstant = fineStructureConstant**2 * electronMass * c / (2 * h)
rydbergConstantHydrogen = rydbergConstant * protonMass / (electronMass + protonMass)
earthRadius = Quantity(6371, km)
sunRadius = solarRadius = Quantity(696342, km)
moonRadius = Quantity(1737.1, km)
earthMass = Quantity(5.97219e24, kg)
sunMass = solarMass = Quantity(1.989e30, kg)
moonMass = Quantity(7.34767309e22, kg)
earthMoonDistance = Quantity(384400, km)
earthSunDistance = Quantity(149597887.5, km)
gravitationalAcceleration = Quantity(9.80665, m / s**2)
zeroCelsius = Quantity(273.15, K)
atmosphericPressure = Quantity(101325, Pa)
waterDensity = Quantity(1000, kg/m**3)

planckLength = (hbar * G / c**3) ** _Fraction(1, 2)
planckMass = (hbar * c / G) ** _Fraction(1, 2)
planckTime = planckLength / c
planckCharge = e / fineStructureConstant**(1/2)
planckEnergy = planckMass * c**2
planckTemperature = planckEnergy / kB
planckArea = planckLength ** 2
planckVolume = planckLength ** 3
planckMomentum = planckMass * c
planckForce = planckEnergy / planckLength
planckPower = planckEnergy / planckTime
planckDensity = planckMass / planckVolume
planckPressure = planckForce / planckArea
planckVoltage = planckEnergy/ planckCharge
planckAcceleration = c / planckTime

minute = Unit(60 * s)
hr = hour = Unit(60 * minute)
d = day = Unit(24 * hour)
w = week = Unit(7 * d)
y = year = Unit(365.2422 * d)
Å = Angstrom = Unit(1e-10 * m)
au = AU = Unit(earthSunDistance)
ly = Unit(9.4605e15 * m)

pc = Unit(3.0857e16 * m)
kpc = Unit(1e3 * pc)
Mpc = Unit(1e6 * pc)
Gpc = Unit(1e9 * pc)

ha = Unit(1e4 * m**2)

b = Unit((1e-14 * m)**2)
mb = Unit(1e-3 * b)
µb = microbarn = Unit(1e-6 * b)
nb = Unit(1e-9 * b)
pb = Unit(1e-12 * b)
fb = Unit(1e-15 * b)

L = Unit(1e-3 * m**3)
dL = Unit(1e-1 * L)
cL = Unit(1e-2 * L)
mL = Unit(1e-3 * L)

g = Unit(1e-3 * kg)
t = Unit(1e3 * kg)
u = Da = Unit(1.660539040e-27 * kg)

eV = Unit(e * V)
meV = Unit(1e-3 * eV)
keV = Unit(1e3 * eV)
MeV = Unit(1e6 * eV)
GeV = Unit(1e9 * eV)
TeV = Unit(1e12 * eV)

atm = Unit(atmosphericPressure)

Ci = Unit(3.7e10 * Bq)
mCi = Unit(1e-3 * Ci)
µCi = microCurie = Unit(1e-6 * Ci)
nCi = Unit(1e-9 * Ci)
pCi = Unit(1e-12 * Ci)

Jy = Unit(1e-26 * W * m**-2 * Hz)
mJy = Unit(1e-3 * Jy)
µJy = microJansky = Unit(1e-6 * Jy)
nJy = Unit(1e-9 * Jy)
pJy = Unit(1e-12 * Jy)
kJy = Unit(1e3 * Jy)

inch = Unit(2.54e-2 * m)
ft = Unit(12 * inch)
yd = Unit(3 * ft)
mi = Unit(1760 * yd)
tsp = Unit(4.92892e-3 * L)
tbsp = Unit(3 * tsp)
floz = Unit(2 * tbsp)
cp = Unit(8 * floz)
pt = Unit(2 * cp)
qt = Unit(2 * pt)
gal = Unit(4 * qt)
oz = Unit(28.3495 * g)
lb = Unit(16 * oz)
cal = Unit(4.19 * J)
kcal = Unit(1e3 * cal)