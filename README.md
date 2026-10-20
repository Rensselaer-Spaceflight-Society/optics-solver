# Geometric Optics Sim

### Note
This *is* a work in progress, so a lot of it is in the works of being implemented.

## Section 1: Components
The structure of setting up a simuation with this library is as follows:

1. Initiate a `Stage` object. This sets any geometric attributes of the space you're working in and contains all the optical elements within separate `Arm`s. The stage object itself is just a container, but may include some sort of metric, if you're working in nonstandard coordinates. The `Arm` objects are moreso where the calculations happen - they contain each individual element and are a handy way of saying:

    >"We have X separate parts of the set, we know they're connected, but we don't want to worry about alignment across the entire setup because that gets messy"

    \- Someone Probably
2. Create a bunch of `Element`s and place them within your setup however you'd like. Each element needs to have a position and rotation, which both have room for N-dimensions in case you need additional coordinates, such as time, angle (e.g. for a pendulum), or even something like polarization.

3. As part of your `Element`s add a form of `DiagnosticElement` somewhere along your beam, and that determines your output.

4. Run via the stage and any relevant graphs will be plotted.

## Section 2: Elements
- `Element`: the generic class of elements serves as a parent for every other element class. This will contain position, rotation, and then also generate a `BoundingBox` which allows checking for overlap between elements.

```python

# Use CoordinateVector class for position
example_position = CoordinateVector(3, 0, 0, 0)
# EulerAngle extends CoordinateVector for rotation
example_rotation = EulerAngle(3, 0, 0, 0)
# create some generic box element
example_element = Element(example_position, example_rotation, 1, 1, 1)

# arg types are :
#   position: CoordinateVector || children
#   rotation: CoordinateVector || children
#   *parameters for whatever
```
Elements will have some function within them that tells the incoming light from a given source how to behave after it passes through. This will be done via a matrix multiplication (where a ray is represented by a 2x2), with gaussian beams it may become a little more complicated.

## Section 3: Sources/Stage/Diagnostics
All of these have been lumped into one section because they essentially serve similar purposes. The source would just be its own class that tells all the elements what's coming in (wavelength, power, coherence time, the works). In order to keep the source as simple as possible, if your source has some kind of initial divergence/aperture you would need to add the class `ApertureElement` to account for that. 

I mentioned the stage working for different coordinates, currently what's in the works is to add a metric function that defines distances and spatial behavior in the system. So for example, could add a class that extends `Stage` called `CurvilinearStage` or `SphericalStage` or even `HyperbolicStage` depending on how we want everything to work. Though a lot of it is generalizing and doesn't need real use.

## Section 4: Experimental Use
Now, simulating is all well and good, but what if we want to (in an actual experiment) find the behavior of some unknown element? I present to you: the `EmpiricalElement` (extends `Element`)! This would theoretically be an element that intakes a dataset based on some form of coordinates from a `CoordinateVector` in the system and some result that would be attached to a `EmpiricalDiagnostic` - which would then be able to calculate the behavior of the element very generally without much effort. What's nice is this can be further simplified if there's already an implemented element that acts *kinda* like what you already have in your system.

### 4.1: Goniometer Example:
So we have the following general setup:
- Light Source
- Sample that rotates
- Sensor

In terms of basic setup we then need:
- CoordinateVectors with 4 coordinates (x, y, z, time).
- An `EmpiricalMirrorElement`, because if we already have an idea of how a mirror behaves (namely fresnel coefficients for reflection), we can more easily get what we need from the data.
- And lastly a `DiagnosticElement` (or one of its children).

Of course, if there are other elements present, such as a fiber optic that essentially acts like an aperture, that can be added with very little extra work.

## Section 5: Extras/Future Stuff
### 5.1: Running Demos:
`python -m demos.one-lens` or more generically: `python -m demos.\[demo_name]` without the `.py`.


### 5.2: Ideas for Future Work
Ideas I have for the future after implementing all the previous plus more:

- Write a demo showing how having n-dimensional things can create a transient system

## OVERALL TODO LIST
- [x] Calculating a hull generator algorithm + basic element class
- [x] Implement nongeneric classes (kinda done)
    - `LensElement`, `MirrorElement`, and `ApertureElement` as a starting point
    - After that look into `DispersiveElement` which may just be a generic thing that causes some attenuation or dispersion
- [x] Implement bounding box Hull class
    - Auxiliary method: find shape overlap between hulls
- [x] Implement Diagnostics and Running and Light Source
- [ ] write a demo for a transient system and make sure it works
- [ ] Write more in-depth docs for each class
