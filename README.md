# Carburizing-in-Python-with-pydiffusion-library
Simple script for multi-step simulation of a carburizing process

This script is a simple way to simulate a multi-step carburizing process. The pydiffusion library is used, a modifified version of sph_sim is included. Carburizing is a well-known and well 
established industrial process. Performed all over the world for applications like gears and shafts, some variations exist. All of these have in common, that elemental carbon is diffused into
the part at rather high temperatures. The standard methods of the pydiffusion library cannot do that, so we applied a slight modification to the sph_sim method.
With a simple modification, just one line, we can turn it into a function, that can simulate the carbon source at the surface. This source potential can easily made time dependent etc.
In the included example, the first carburizing step is followed by a diffusion step with the standard method sph_sim. A second carburizing and diffusion step follow.
In principle, this setup can be used for any diffusing species, steps can be added/removed to your liking.
The diffusion coefficient in the example is a f(T). 
It thus follows, that possibilities with this script are nearly limitless.

Output of first plot:

![Carburizing+Diffusion_steps1](https://github.com/emefff/Carburizing-in-Python-with-pydiffusion-library/assets/89903493/f0539904-af4b-4a05-b93c-88226c9d82a3)

Output of second plot:

![Carburizing+Diffusion_steps2](https://github.com/emefff/Carburizing-in-Python-with-pydiffusion-library/assets/89903493/590baf6e-fcc8-4e34-99d8-8ba50178cfad)
