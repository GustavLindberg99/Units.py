# -*- coding: utf-8 -*-
# units.py version 2.0.1 by Gustav Lindberg
# https://github.com/GustavLindberg99/Units.py
# https://pypi.org/project/unitspy/

from __future__ import annotations
from fractions import Fraction
from math import isclose as _real_isclose
from math import pi
from typeguard import typechecked
import re as _re

def _isclose(a: complex, b: complex) -> bool:
    return _real_isclose(a.real, b.real) and _real_isclose(a.imag, b.imag)

#Otherwise it will cause errors when using Quantity as type hints in Unit
class Quantity:
    pass

#Don't add @typechecked to stuff that doesn't need it, otherwise it's bad for performance

class Unit:
    _kilograms: int | Fraction
    _meters: int | Fraction
    _seconds: int | Fraction
    _amperes: int | Fraction
    _kelvins: int | Fraction
    _moles: int | Fraction
    _multiple: float
    
    @typechecked
    def __init__(this, other: Unit | Quantity | float | None = None):
        if isinstance(other, Unit):
            this._kilograms = other._kilograms
            this._meters = other._meters
            this._seconds = other._seconds
            this._amperes = other._amperes
            this._kelvins = other._kelvins
            this._moles = other._moles
            m = other._multiple
            this._multiple = m
        elif isinstance(other, Quantity):
            this.__init__(other._unit)
            this._multiple *= other._number
        else:
            this._kilograms = 0
            this._meters = 0
            this._seconds = 0
            this._amperes = 0
            this._kelvins = 0
            this._moles = 0
            if other == None:
                this._multiple = 1
            else:
                this._multiple = other
        
    def __copy__(this):
        return Unit(this)
    
    def __deepcopy__(this):
        return Unit(this)
    
    
    @typechecked
    def __imul__(this, other: Unit | Quantity | float) -> Unit:
        if isinstance(other, Quantity):
            other = Unit(other)
        if isinstance(other, Unit):
            this._kilograms += other._kilograms
            this._meters += other._meters
            this._seconds += other._seconds
            this._amperes += other._amperes
            this._kelvins += other._kelvins
            this._moles += other._moles
            this._multiple *= other._multiple
        else:
            this._multiple *= other
        return this
    
    def __mul__(this, other: Unit | Quantity | complex) -> Unit | Quantity:
        result: Quantity = Quantity(this)
        result *= other
        if isinstance(other, Unit):
            return Unit(result)
        elif isinstance(other, Quantity):
            return result.toUnit(this * other._unit)
        else:
            return result.toUnit(this)
    
    def __rmul__(this, other: Unit | Quantity | complex) -> Unit | Quantity:
        return this * other
    
    
    @typechecked
    def __ipow__(this, other: float | Fraction) -> Unit:
        other = Fraction(other)
        this._kilograms *= other
        this._meters *= other
        this._seconds *= other
        this._amperes *= other
        this._kelvins *= other
        this._moles *= other
    
        this._multiple = float(this._multiple ** other)
        return this
        
    def __pow__(this, other: float | Fraction) -> Unit:
        result: Unit = Unit(this)
        result **= other
        return result
    
    
    def __idiv__(this, other: Unit | Quantity | float) -> Unit:
        this *= other ** -1
        return this
        
    def __truediv__(this, other: Unit | Quantity | complex) -> Unit | Quantity:
        return this * other ** -1
    
    def __rtruediv__(this, other: Unit | Quantity | complex) -> Unit | Quantity:
        return other * this ** -1
    
    
    @typechecked
    def isMultiple(this, other: Unit) -> bool:
        return (
            this._kilograms == other._kilograms and
            this._meters == other._meters and
            this._seconds == other._seconds and
            this._amperes == other._amperes and
            this._kelvins == other._kelvins and
            this._moles == other._moles
        )
        
    def isSiUnit(this) -> bool:
        return _isclose(this._multiple, 1.0)
    
    def isPlanckUnit(this) -> bool:
        return Quantity(1, this).toPlanckUnits()._unit == this
            
    
    @typechecked
    def __eq__(this, other: Unit | Quantity | complex) -> bool:
        if not isinstance(other, Unit):
            return this == Unit(other)
        return this.isMultiple(other) and _isclose(this._multiple, other._multiple)
        
    def __ne__(this, other: Unit | Quantity | complex) -> bool:
        return not(this == other)
    
    def __nonzero__(this) -> bool:
        return this._multiple != 0
    
    
    def __neg__(this) -> Unit:
        return this * -1
    
    def __pos__(this) -> Unit:
        return this
    
    
    def __hash__(this) -> int:
        return hash((
            this._kilograms,
            this._meters,
            this._seconds,
            this._amperes,
            this._kelvins,
            this._moles,
            this._multiple
        ))
        
    def __int__(this) -> int:
        return int(float(this))
        
    def __float__(this) -> float:
        if not(this.isMultiple(dimensionless)):
            raise ValueError("Only units that are multiples of dimensionless can be converted to float or complex, not " + str(this))
        return float(this._multiple)
    
    def __repr__(this) -> str:
        if this == dimensionless:
            return "Unit(1)"
        elif this._moles == 1 and this != mol:
            thisTimesNA: str = str(Unit(this * NA))
            if(not thisTimesNA[0].isdigit() or thisTimesNA == "1"):
                return thisTimesNA
            return str(this / mol) + "*mol"
        elif this._moles == -1 and this != mol:
            thisDividedByNA: str = str(Unit(this / NA))
            if not thisDividedByNA[0].isdigit() or thisDividedByNA == "1":
                return thisDividedByNA
            result: str = str(this * mol)
            if('/' in result):
                result += "/mol"
            else:
                result += "*mol**-1"
            return result
        elif(this == t):
            return "t"
        elif(this == minute):
            return "min"
        elif(this == minute**-1):
            return "min**-1"
        elif(this == hour):
            return "hr"
        elif(this == hour**-1):
            return "hr**-1"
        elif(this == d):
            return "d"
        elif(this == d**-1):
            return "d**-1"
        elif(this == w):
            return "w"
        elif(this == w**-1):
            return "w**-1"
        elif(this == y):
            return "y"
        elif(this == y**-1):
            return "y**-1"
        elif(this == Å):
            return "Å"
        elif(this == Å**2):
            return "Å**2"
        elif(this == Å**3):
            return "Å**3"
        elif(this == au):
            return "au"
        elif(this == au**2):
            return "au**2"
        elif(this == au**3):
            return "au**3"
        elif(this == ly):
            return "ly"
        elif(this == ly**2):
            return "ly**2"
        elif(this == ly**3):
            return "ly**3"
        elif(this == ha):
            return "ha"
        elif(this == b):
            return "b"
        elif(this == inch):
            return "in"
        elif(this == inch**2):
            return "in**2"
        elif(this == inch**3):
            return "in**3"
        elif(this == ft):
            return "ft"
        elif(this == ft**2):
            return "ft**2"
        elif(this == ft**3):
            return "ft**3"
        elif(this == yd):
            return "yd"
        elif(this == yd**2):
            return "yd**2"
        elif(this == yd**3):
            return "yd**3"
        elif(this == mi):
            return "mi"
        elif(this == mi**2):
            return "mi**2"
        elif(this == mi**3):
            return "mi**3"
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
            return "deg"
        elif(this == Unit(e)):
            return "Unit(e)"
        elif(this == Unit(G)):
            return "Unit(G)"
        elif(this == Unit(hbar)):
            return "Unit(hbar)"
        elif(this == Unit(c)):
            return "Unit(c)"
        elif(this == Unit(epsilon0)):
            return "Unit(epsilon0)"
        elif(this == Unit(mu0)):
            return "Unit(mu0)"
        elif(this == Unit(planckEnergy)):
            return "Unit(E_p)"
        unitToString: dict[Unit, str] = dict()
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
        unitToString[Ohm] = "Ohm"
        unitToString[S] = "S"
        unitToString[Wb] = "Wb"
        unitToString[T] = "T"
        unitToString[H] = "H"
        basicUnits: dict[Unit, str] = dict(unitToString)
        unitToString[L] = "L"
        unitToString[g] = "g"
        unitToString[eV] = "eV"
        unitToString[Unit(eV / c)] = "eV/c"
        unitToString[Unit(eV / c**2)] = "eV/c**2"
        unitToString[mol] = "mol"
        unitToString[Da] = "Da"
        unitToString[cal] = "cal"
        unitToString[pc] = "pc"
        for unit in list(basicUnits.keys()) + [g, L, Unit(eV), Unit(eV / c), Unit(eV / c**2), mol, Da, cal, pc]:
            unitToString[Unit(1e-1 * unit)] = "d" + unitToString[unit]
            unitToString[Unit(1e-2 * unit)] = "c" + unitToString[unit]
            unitToString[Unit(1e-3 * unit)] = "m" + unitToString[unit]
            unitToString[Unit(1e-6 * unit)] = "micro" + unitToString[unit]
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
        
        def powerString(unit: str, power: int | Fraction) -> str:
            if power == 1:
                return "{}*".format(unit)
            elif power == 0:
                return ""
            elif isinstance(power, Fraction) and power.denominator != 1:
                return "{}**({})*".format(unit, power)
            return "{}**{}*".format(unit, power)
        
        try:
            if this.isPlanckUnit():
                result: str = powerString("m_p", this._kilograms) + powerString("l_p", this._meters) + powerString("t_p", this._seconds - this._amperes) + powerString("q_p", this._amperes) + powerString("T_p", this._kelvins)
                if this._moles != 0 and result == "":
                    return "Unit(1)"
                return "Unit({})".format(result[:-1])
        except OverflowError:
            pass
            
        if this in unitToString.keys():
            return unitToString[this]
        for unit in basicUnits.keys():
            if this.isMultiple(unit):
                return "Unit({}*{})".format(this._multiple, unitToString[unit])
        if this == dimensionless:
            return "Unit(1)"
        if this.isMultiple(dimensionless):
            return "Unit({})".format(this._multiple)
        for unit in unitToString.keys():
            for i in range(-3, 4):
                if this == unit ** i:
                    return unitToString[unit] + "**" + str(i)
        for unit in basicUnits.keys():
            for i in range(-3, 4):
                if(this.isMultiple(unit ** i)):
                    return str(this._multiple) + "*" + unitToString[unit] + "**" + str(i)
        units = [m, s, J, eV, V, K, F, H, W, C]
        for unit1 in units:
            for unit2 in units:
                for i in [1, 2, 3]:
                    if(this == unit1 * unit2**i):
                        return unitToString[unit1] + "*" + unitToString[unit2] + ("" if i == 1 else "**" + str(i))
                    elif(this == unit1 / unit2**i):
                        return unitToString[unit1] + "/" + unitToString[unit2] + ("" if i == 1 else "**" + str(i))
        
        result: str = (("" if this.isSiUnit() else "{}*".format(this._multiple))
            + powerString("kg", this._kilograms)
            + powerString("m", this._meters)
            + powerString("s", this._seconds)
            + powerString("A", this._amperes)
            + powerString("K", this._kelvins)
            + powerString("mol", this._moles)
        )[:-1]
        if not this.isSiUnit():
            return "Unit({})".format(result)
        return result
    
    def __str__(this) -> str:
        def rawSuperscriptToUnicode(rawSuperscript: _re.Match) -> str:
            result: str = ""
            for character in rawSuperscript.group(0):
                match character:
                    case '-':
                        result += "\u207b"
                    case '0':
                        result += "\u2070"
                    case '1':
                        result += "\u00b9"
                    case '2':
                        result += "\u00b2"
                    case '3':
                        result += "\u00b3"
                    case '4':
                        result += "\u2074"
                    case '5':
                        result += "\u2075"
                    case '6':
                        result += "\u2076"
                    case '7':
                        result += "\u2077"
                    case '8':
                        result += "\u2078"
                    case '9':
                        result += "\u2079"
            return result
        
        result: str = this.__repr__()
        result = _re.sub(r"(Unit|Quantity)\((.+)\)$", lambda match: match.group(2), result)
        result = _re.sub(r"\*\*(\(-?[0-9]+/[0-9]+\)|-?[0-9]+)", rawSuperscriptToUnicode, result)
        result = _re.sub(r"(?<=[0-9])\*", "", result)
        result = result.replace("deg", "\u00b0")
        result = result.replace("hbar", "\u0127")
        result = result.replace("epsilon0", "\u03b5\u2080")
        result = result.replace("mu0", "\u03bc\u2080")
        result = result.replace("Ohm", "\u03a9")
        result = result.replace("_p", "\u209a")
        result = result.replace("micro", "\u03bc")
        result = result.replace("*", "·")
        return result

#Inherit from the dummy class defined above so that it doesn't cause errors of the type "object of type Quantity is not an instance of Quantity"
class Quantity(Quantity):
    _unit: Unit
    _number: complex
    
    @typechecked
    def __init__(this, number: complex | Quantity | Unit = 1.0, unit: Unit | Quantity = Unit()):
        this._unit = Unit(unit)
        if isinstance(number, Quantity):
            this._unit *= number._unit
            this._number = number._number
        elif isinstance(number, Unit):
            this._unit = Unit(number)
            this._number = 1
        else:
            this._number = number
    
    
    def __copy__(this):
        result = Quantity()
        result._number = this._number
        result._unit = this._unit
        return result
    
    def __deepcopy__(this):
        return Quantity(this)
        
        
    @typechecked
    def __imul__(this, other: Quantity | Unit | complex) -> Quantity:
        if isinstance(other, Quantity):
            this._number *= other._number
            this._unit *= other._unit
        elif isinstance(other, Unit):
            this._unit *= other
        else:
            this._number *= other
        return this
    
    def __mul__(this, other: Quantity | Unit | complex) -> Quantity | complex:
        result = Quantity(this)
        result *= other
        return result.dimensionlessToFloat()
    
    def __rmul__(this, other: Quantity | Unit | complex) -> Quantity | complex:
        return this * other
    
    
    @typechecked
    def __ipow__(this, other: complex | Fraction) -> Quantity:
        this._number **= other
        if isinstance(this._number, Fraction):
            this._number = float(this._number)
        this._unit **= other
        return this
        
    def __pow__(this, other: complex | Fraction) -> Quantity | complex:
        result: Quantity = Quantity(this)
        result **= other
        return result.dimensionlessToFloat()
    
    
    def __idiv__(this, other: Quantity | Unit | complex) -> Quantity:
        this *= other ** -1
        return this
    
    def __truediv__(this, other: Quantity | Unit | complex) -> Quantity | complex:
        return this * other ** -1
    
    def __rtruediv__(this, other: Quantity | Unit | complex) -> Quantity | complex:
        return other * this ** -1
    
    
    def __neg__(this) -> Quantity | complex:
        return (this * -1).dimensionlessToFloat()
    
    def __pos__(this) -> Quantity | complex:
        return this.dimensionlessToFloat()
    
    
    @typechecked
    def __iadd__(this, other: Quantity | complex) -> Quantity:
        if other == 0 and Quantity(other)._unit.isMultiple(dimensionless):
            return this
        if this == 0 and Quantity(this)._unit.isMultiple(dimensionless):
            this._number = other._number
            this._unit = other._unit
            return this
        other = Quantity(other)
        if not this._unit.isMultiple(other._unit):
            raise ValueError("Cannot add two quantities with different units: {} and {}".format(this._unit, other._unit))
        this._number += other._number * Quantity(other._unit / this._unit).dimensionlessToFloat()
        return this
        
    def __add__(this, other: Quantity | complex) -> Quantity | complex:
        result: Quantity = Quantity(this)
        result += other
        return result.dimensionlessToFloat()
    
    def __radd__(this, other: Quantity | complex) -> Quantity | complex:
        return this + other
    
    
    def __isub__(this, other: Quantity | complex) -> Quantity:
        this += (-other)
        return this
    
    def __sub__(this, other: Quantity | complex) -> Quantity | complex:
        return this + (-other)
    
    def __rsub__(this, other: Quantity | complex) -> Quantity | complex:
        return -(this - other)
    
    
    @typechecked
    def __eq__(this, other: Quantity | Unit | complex) -> bool:
        if other == 0:
            return this._number == 0 or this._unit._multiple == 0
        elif isinstance(other, Quantity):
            return this._unit.isMultiple(other._unit) and _isclose((this._number * this._unit._multiple) / (other._number * other._unit._multiple), 1)
        elif isinstance(other, Unit):
            return this == Quantity(1, other)
        elif this._unit.isMultiple(dimensionless):
            return complex(this) == other
        else:
            return False
        
    def __neq__(this, other: Quantity | Unit | complex) -> bool:
        return not(this == other)
    
    def __nonzero__(this) -> bool:
        return this != 0
    
    
    @typechecked
    def __lt__(this, other: Quantity | complex) -> bool:
        if this == 0 and Quantity(this)._unit.isMultiple(dimensionless):
            return 0 <= other._number
        if other == 0 and Quantity(other)._unit.isMultiple(dimensionless):
            return this._number <= 0
        other = Quantity(other)
        if not this._unit.isMultiple(other._unit):
            raise ValueError("Cannot compare two quantities with different units: {} and {}".format(this._unit, other._unit))
        return this._number * this._unit._multiple < other._number * other._unit._multiple
    
    def __le__(this, other: Quantity | complex) -> bool:
        return this < other or this == other
    
    def __gt__(this, other: Quantity | complex) -> bool:
        return not(this <= other)
    
    def __ge__(this, other: Quantity | complex) -> bool:
        return not(this < other)
    
    
    def __hash__(this) -> int:
        return hash((this._number, this._unit))
    
    def __int__(this) -> int:
        return int(float(this))
        
    def __float__(this) -> float:
        if complex(this).imag != 0:
            raise ValueError("Cannot convert complex Quantity to float or int")
        return complex(this).real
    
    def __complex__(this) -> complex:
        return complex(this._number) * float(this._unit)
    
    
    def __repr__(this) -> str:
        if this._unit == dimensionless:
            return "Quantity({})".format(this._number)
        elif not this._unit.__repr__()[0].isdigit() and this._unit.__repr__()[0] != '.' and this._unit.__repr__()[0] != '-':
            return this._number.__repr__() + "*" + _re.sub("^[0-9.-]+", "", this._unit.__repr__())
        elif this._unit.isMultiple(dimensionless):
            return "Quantity({})".format(this._number * float(this._unit))
        else:
            return (this._number * this._unit._multiple).__repr__() + "*" + _re.sub("^[0-9.-]+", "", Unit(this._unit / this._unit._multiple).__repr__())
    
    def __str__(this) -> str:
        return Unit.__str__(this)
        
    
    @typechecked
    def toUnit(this, unit: Unit) -> Quantity:
        if(this._unit._kilograms != unit._kilograms or this._unit._meters != unit._meters or this._unit._seconds != unit._seconds or this._unit._amperes != unit._amperes or this._unit._kelvins != unit._kelvins):
            raise ValueError("Cannot convert quantity in unit {} to quantity in unit {}".format(this._unit, unit))
        result: Quantity = Quantity(this._number * this._unit._multiple / unit._multiple, unit)
        if(this._unit._moles != unit._moles):
            result *= (NA * mol) ** (this._unit._moles - unit._moles)
        return result
    
    def toSiUnits(this) -> Quantity:
        unit: Unit = Unit(this._unit)
        unit._multiple = 1
        return this.toUnit(unit)
    
    def toPlanckUnits(this) -> Quantity:
        unit: Unit = Unit(planckMass ** this._unit._kilograms * planckLength ** this._unit._meters * planckTime ** this._unit._seconds * (planckCharge / planckTime) ** this._unit._amperes * planckTemperature ** this._unit._kelvins)
        return this.toUnit(unit)
    
    def molesToDimensionless(this) -> Quantity:
        unit: Unit = Unit(this._unit)
        unit._moles = 0
        return this.toUnit(unit)
    
    def isDimensionless(this) -> bool:
        return this._unit.isMultiple(dimensionless)
    
    def dimensionlessToFloat(this) -> Quantity | complex:
        if(this.isDimensionless()):
            if this._number.imag == 0:
                return float(this)
            else:
                return complex(this)
        else:
            return this
    
    def __abs__(this) -> Quantity:
        return Quantity(abs(this._number), this._unit)
    
    def __round__(this) -> Quantity:
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
deg = degree = degrees = Unit(pi / 180)
arcmin = Unit(degree / 60)
arcsec = Unit(degree / 3600)
radian = dimensionless

cm = Unit(1e-2 * m)
mm = Unit(1e-3 * m)
µm = micrometer = microm = Unit(1e-6 * m)
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
µA = microAmpere = microA = Unit(1e-6 * A)
nA = Unit(1e-9 * A)
pA = Unit(1e-12 * A)
fA = Unit(1e-15 * A)

mF = Unit(1e-3 * F)
µF = microFarad = microF = Unit(1e-6 * F)
nF = Unit(1e-9 * F)
pF = Unit(1e-12 * F)

mOhm = milliOhm = Unit(1e-3 * Ohm)
kOhm = kiloOhm = Unit(1e3 * Ohm)
MOhm = megaOhm = Unit(1e6 * Ohm)

ms = Unit(1e-3 * s)
µs = microsecond = micros = Unit(1e-6 * s)
ns = Unit(1e-9 * s)
ps = Unit(1e-12 * s)
fs = Unit(1e-15 * s)

mg = Unit(1e-6 * kg)
µg = microgram = microg = Unit(1e-9 * kg)
ng = Unit(1e-12 * kg)

kBq = Unit(1e3 * Bq)
MBq = Unit(1e6 * Bq)
kHz = Unit(1e3 * Hz)
MHz = Unit(1e6 * Hz)
GHz = Unit(1e9 * Hz)
THz = Unit(1e12 * Hz)

mC = Unit(1e-3 * C)
µC = microCoulomb = microC = Unit(1e-6 * C)
nC = Unit(1e-9 * C)
pC = Unit(1e-12 * C)
fC = Unit(1e-15 * C)

i = 1j
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

planckLength = l_p = (hbar * G / c**3) ** Fraction(1, 2)
planckMass = m_p = (hbar * c / G) ** Fraction(1, 2)
planckTime = t_p = planckLength / c
planckCharge = q_p = e / fineStructureConstant**(1/2)
planckEnergy = E_p = planckMass * c**2
planckTemperature = T_p = planckEnergy / kB
planckArea = A_p = planckLength ** 2
planckVolume = V_p = planckLength ** 3
planckMomentum = p_p = planckMass * c
planckForce = F_p = planckEnergy / planckLength
planckPower = planckEnergy / planckTime
planckDensity = rho_p = planckMass / planckVolume
planckPressure = planckForce / planckArea
planckVoltage = planckEnergy/ planckCharge
planckAcceleration = a_p = c / planckTime

minute = Unit(60 * s)
hr = hour = Unit(60 * minute)
d = day = Unit(24 * hour)
w = week = Unit(7 * d)
y = year = yr = Unit(365.2422 * d)
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
µb = microbarn = microb = Unit(1e-6 * b)
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
µCi = microCurie = microCi = Unit(1e-6 * Ci)
nCi = Unit(1e-9 * Ci)
pCi = Unit(1e-12 * Ci)

Jy = Unit(1e-26 * W * m**-2 * Hz)
mJy = Unit(1e-3 * Jy)
µJy = microJansky = microJy = Unit(1e-6 * Jy)
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
