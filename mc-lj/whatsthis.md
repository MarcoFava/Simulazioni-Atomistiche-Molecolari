# Lab task:
_The final goals of next two lab sessions are to (i) develop a minimal but flexibile Monte Carlo code to simulate classical many-body systems and (ii) use it compute the equation of state of a Lennard-Jones fluid_

There will be a **system** *composed* by N **particles** in 3 **dimensions**, confined in a **box** of given **sides** with **perodic boundary conditions** and *effected by* a **potential** (for us it'll be L-J pot).  
Each particle has a **position** *starting* from an **initial configuration** (crystal).  
The system has a certain **temperature**, **density**, **pressure**, ...  

We will *evolve* the system by *using* **montecarlo methods**...

---
- **OBJECTS and ATTRIBUTES**
    - system
    - particles
    - dimensions
    - box
    - sides
    - pbc
    - potential
    - position
    - initial configuration
    - temperature
    - density
    - pressure

- **METHODS and OBJECT RELATIONSHIPS**
    - system **composed by**
    - **confined** in a box
    - **effected by** a potential
    - **starting from** an initial configuration
    - **evolve** the system
    - **using** montecarlo methods