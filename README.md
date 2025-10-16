# SPH Water Discharge Simulation
This project implements a 2D water discharge simulation using Smoothed Particle Hydrodynamics (SPH). The simulation models a tank of water under gravitational force with a dynamically opening hole, allowing visualization of fluid dynamics in confined spaces.

**Key Features**
  
  SPH Method: Uses particle-based approach to simulate incompressible fluid behavior
  
  Dynamic Scenario: Initially simulates water settling in a closed tank, then opens a hole in the right wall to observe discharge

  Physics Modeling:
  
    *Implements Navier-Stokes equations for fluid dynamics
    
    *Applies gravitational forces
    
    *Handles boundary interactions with walls
    
    *Computes particle density with kernel smoothing functions
    
    *Uses equation of state (EOS) for pressure calculation
    
    *Incorporates XSPH velocity correction for stability

**Technical Implementation**

  *Efficient nearest-neighbor search using grid-based linked list approach
  
  *Leap-frog integration scheme for time evolution
  
  *Custom boundary handling to prevent particle penetration
  
  *Adaptive density normalization for improved accuracy
  
  *C implementation with optimized data structures

**Visualization**

The simulation outputs time-series data files that can be visualized to observe the evolution of the water discharge process, from initial stabilization under gravity to the dynamic flow through the opening.

This project demonstrates fundamental concepts in computational fluid dynamics using particle methods, with applications in engineering simulations and computer graphics.
