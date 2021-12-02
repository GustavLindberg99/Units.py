# units.py
Units.py is a Python library that makes it easy to work with units in physics without having to worry about unit conversions. For example, if you have an energy in electronvolts and you have the value of Planck's constant in Joules times seconds, normally you would have to convert one of those values to match the units of the other one. But with units.py, you don't have to do that. You can easily multiply two values in different units and get the correct result.

# Setup
To start using units.py, download [this file](https://raw.githubusercontent.com/Gustav-Lindberg/Units.py/main/units.py) and save it in the same folder as your Python code. Then add this line to the beginning of your Python code:

```
import units
```

If you want it to be available for all projects, you can also save the file in the folder containing Python libraries. On Windows with Anaconda 3 this is `C:\ProgramData\Anaconda3\Lib`, and on Linux with Python 3.8 this is `/usr/lib/python3.8`.

# Usage
## Introduction
Units.py creates a `Quantity` class. A `Quantity` object represents a physical quantitiy, for example a length, a time or a mass. To create a `Quantity` object, simply multiply a number by a unit provided by units.py. For example, to create a quantity representing 3 meters, do this:

```
threeMeters = 3 * units.m
```

You can multiply and divide `Quantity` objects just like numbers, and you can add, subtract and compare two `Quantity` objects as long as they have compatible units (for example seconds and days are compatible units as they're both units of time, seconds and meters are not as they measure different things). Example:

```
myLength = 10 * units.m
myTime = 5 * units.s
print(myLength / myTime)    #Result: 2m/s
```

## Units
Most units are available, and their names are often their SI abbreviations. Here is a list of the most common units:

Units.py name (should be preceded by `units.`) | Unit
--- | ---
`g` | [Gram](https://en.wikipedia.org/wiki/Gram)
`kg` | [Kilogram](https://en.wikipedia.org/wiki/Kilogram)
`m` | [Meter](https://en.wikipedia.org/wiki/Metre)
`Å` or `Angstrom` (the two are equivalent) | [Ångström](https://en.wikipedia.org/wiki/Angstrom)
`au` | [Astronomical unit](https://en.wikipedia.org/wiki/Astronomical_unit)
`ly` | [Light year](https://en.wikipedia.org/wiki/Light-year)
`pc` | [Parsec](https://en.wikipedia.org/wiki/Parsec)
`s` | [Second](https://en.wikipedia.org/wiki/Second)
`minute` | [Minute](https://en.wikipedia.org/wiki/Minute)
`hr` | [Hour](https://en.wikipedia.org/wiki/Hour)
`d` | [Day](https://en.wikipedia.org/wiki/Day)
`y` | [Year](https://en.wikipedia.org/wiki/Year)
`A` | [Ampere](https://en.wikipedia.org/wiki/Ampere)
`K` | [Kelvin](https://en.wikipedia.org/wiki/Kelvin)
`Hz` | [Hertz](https://en.wikipedia.org/wiki/Hertz)
`N` | [Newton](https://en.wikipedia.org/wiki/Newton_(unit))
`Pa` | [Pascal](https://en.wikipedia.org/wiki/Pascal_(unit))
`J` | [Joule](https://en.wikipedia.org/wiki/Joule)
`eV` | [Electronvolt](https://en.wikipedia.org/wiki/Electronvolt)
`W` | [Watt](https://en.wikipedia.org/wiki/Watt)
`C` | [Coulomb](https://en.wikipedia.org/wiki/Coulomb)
`V` | [Volt](https://en.wikipedia.org/wiki/Volt)
`F` | [Farad](https://en.wikipedia.org/wiki/Farad)
`Ohm` | [Ohm](https://en.wikipedia.org/wiki/Ohm)
`T` | [Tesla](https://en.wikipedia.org/wiki/Tesla_(unit))

You can also add SI prefixes before some of these units to get multiples of them. For example, write `cm` for centimeter or `ms` for millisecond. For units starting with "micro", you can either precede the unit with `µ` (for example `µm`), or write out the the name of the entire unit (for example `micrometer`).

You can also multiply, divide and take the power of units. For example, you can write `m / s` for meters per second or `m ** 2` for square meters.

### Angles and dimensionless quantities
Units.py also provides a `dimensionless` unit in order to represent [dimensionless quantities](https://en.wikipedia.org/wiki/Dimensionless_quantity). A `Quantity` object with the `dimensionless` unit is meant to work exactly like an object of Python's built-in `float` type.

For angles, there is also a `radians` unit, which is simply an alias of `dimensionless`, and a `degrees` unit, where `1 * units.degrees` is equal to `pi / 180 * units.radians`. That way, you can easily specify angles in degrees, and for example pass that as an argument to trigonometric functions like `numpy.sin`, for example

```
import units
import numpy
print(numpy.sin(90 * units.degree))    #Result: 1.0
```

There is also `units.arcmin` and `units.arcsec`, representing arcminutes and arcseconds respectively.

### Temperatures
Units.py stores temperatures in Kelvins, and can't store temperatures in Celsius. However, it defines `Quantity` object named `zeroCelsius` equal to the value of 0 Celsius in Kelvin, which can be used to convert between Kelvin and Celsius. If you have a temperature in Celsius, add `zeroCelsius` to get the temperature in Kelvin, and if you have a temperature in Kelvin, subtract `zeroCelsius` to get the temperature in Celsius. Example:

```
temperatureInCelsius = 15 * units.K    #Define this quantity in Kelvin even though it really represents a unit in Celsius
temperatureInKelvin = temperatureInCelsius + units.zeroCelsius
print(temperatureInKelvin)    #Result: 288.15K
```

## Physical quantities
Units.py also provides pre-defined `Quantity` objects corresponding to physical constants. These are some common ones:

Units.py name (should be preceded by `units.`) | Quantity
--- | ---
`kB` | [Boltzmann constant](https://en.wikipedia.org/wiki/Boltzmann_constant)
`e` | [Elementary charge](https://en.wikipedia.org/wiki/Elementary_charge)
`G` | [Gravitational constant](https://en.wikipedia.org/wiki/Gravitational_constant)
`h` | [Planck constant](https://en.wikipedia.org/wiki/Planck_constant)
`hbar` | [Reduced Planck constant](https://en.wikipedia.org/wiki/Reduced_Planck_constant)
`c` | [Speed of light](https://en.wikipedia.org/wiki/Speed_of_light)
`epsilon0` | [Vacuum permittivity](https://en.wikipedia.org/wiki/Vacuum_permittivity)
`mu0` | [Vacuum permeability](https://en.wikipedia.org/wiki/Vacuum_permeability)
`me` | [Electron mass](https://en.wikipedia.org/wiki/Electron_rest_mass)
`mp` | [Proton](https://en.wikipedia.org/wiki/Proton) mass
`mn` | [Neutron](https://en.wikipedia.org/wiki/Neutron) mass
`solarMass` | [Solar mass](https://en.wikipedia.org/wiki/Solar_mass)
`pi` | [Pi](https://en.wikipedia.org/wiki/Pi)

## Conversions between units
If you want to output the value of a `Quantity` object in a specific unit, you will need the `toUnit` method. This method takes one argument, which is the unit in which you want to output the quantity. For example, here is how you output the value of one light year in meters:

```
myDistance = 1 * units.ly
print(myDistance.toUnit(units.m))    #Result: 9460500000000000.0m
```

`Quantity` objects also have a `toSiUnits` method which converts it to [SI units](https://en.wikipedia.org/wiki/International_System_of_Units) and a `toPlanckUnits` method which converts it to [Planck units](https://en.wikipedia.org/wiki/Planck_units). These methods take no arguments. Example:

```
myDistance = 1 * units.ly
print(myDistance.toSiUnits())    #Result: 9460500000000000.0m
print(myDistance.toPlanckUnits())    #Result: 5.853346072890479e+50lₚ
```
