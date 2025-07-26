# Navier-2D

A 2D incompressible Navier-Stokes fluid dynamics solver implemented in C++.

---

## Overview

**Navier-2D** is a C++ program that solves the two-dimensional incompressible Navier-Stokes equations to simulate fluid flow. It computes velocity and pressure fields on a 2D grid using finite difference methods.

This project is ideal for learning computational fluid dynamics and experimenting with fluid flow simulations.

---

## Features

- Solves 2D incompressible Navier-Stokes equations  
- Configurable viscosity, grid size, and time steps  
- Supports common boundary conditions  
- Outputs velocity and pressure data for visualization  
- Modular code structure for easy extensions  

---

## Building the Project

Make sure you have a C++ compiler installed (e.g., g++, clang++).

To build:

```bash
g++ -o navier2d navier2d.cpp
