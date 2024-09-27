# Gerchberg and Saxton Algorithm for Phase Retrieval and Projection on Spatial Light Modulator

This repository implements the **Gerchberg-Saxton Algorithm** to retrieve the phase of a desired intensity pattern and project it onto a Spatial Light Modulator (SLM). It allows the user to input beam parameters, generate target images (e.g., multiple spots), and project the calculated *phase mask* or *digital hologram* onto the SLM for beam shaping. This digital hologram is particularly useful in optical setups where precise phase control is required to shape light beams. Here, the digital hologram generates an array of *optical traps*. The file contains, `Multiple_Spto.m`, which generates the desired input pattern for the phase mask. The phase mask is generated using the Gerchberg Saxton Algorithm `GS_ALgo.m`, it is also called a digital hologram.  
This computer-generated phase mask is projected onto the spatial light modulator and the desired pattern is obtained at the focal plane of the objective resulting in the series of *optical traps*. 

**Note:** While using this code, make sure that all `.m` files are together.

## Features
- **Gerchberg-Saxton Phase Retrieval Algorithm**: Iterative algorithm to retrieve the phase information from an input intensity pattern.
- **SLM Projection**: Projects the phase mask onto an SLM for holographic applications.
- **Flexible Input Beam Definition**: Customizable input beam profile and diffracted optical element.
- **Target Image Generation**: Ability to create target intensity patterns, such as multiple spots, based on user input.
- **Random Phase Initialization**: Random phase is applied initially for phase retrieval.

## How it Works

1. **Input Beam Definition**: 
   - The input beam is circular with a defined radius.
   - It is represented as an amplitude field.
  
2. **Gerchberg-Saxton Algorithm**: 
   - The phase of the input beam is initialized randomly.
   - The algorithm iteratively computes the phase mask by forward propagating the wave through Fourier space and adjusting the phase in each iteration to match the desired intensity pattern.
   
3. **Phase Mask Generation**: 
   - The phase mask generated from the algorithm is displayed and can be used for projection onto an SLM.

4. **Projection onto SLM**: 
   - The generated phase mask is projected onto the SLM. The resolution of the SLM must match the size of the generated phase mask.
   - The phase mask is displayed full-screen on the secondary monitor to simulate projection onto the SLM.

## Usage

1. **Input Beam**: 
   - Modify the radius and position of the input beam as needed by editing the code in the input beam section.

2. **Target Image**:
   - You will be prompted to enter the number of spots and their separation when running the code.
   - This creates a custom intensity pattern for the algorithm to shape.

3. **Iterations**:
   - Enter the desired number of iterations for the Gerchberg-Saxton algorithm. Typically, more iterations yield a more accurate phase mask.

4. **Projection**:
   - Ensure that the SLM resolution matches the size of the generated phase mask (typically set to match the resolution of the secondary monitor).
   - The phase mask will be displayed on the secondary monitor (SLM screen).

## Instructions

1. Clone this repository:
   ```bash
   git clone https://github.com/hisaylama/gerchberg-saxton-slm.git
   ```

2. Open the MATLAB script in your MATLAB environment.

3. Run the script. You will be prompted for input values such as the number of spots, separation, focal length, and the number of iterations.

4. The phase mask will be generated and displayed. If an SLM is connected, the phase mask will be projected onto the SLM.

5. The result can be visualized on your SLM or secondary monitor.

## Requirements

- **MATLAB**: The script is designed to run in a MATLAB environment.
- **SLM**: A spatial light modulator with resolution matching the input mask is required for projection.
- **Monitor Setup**: The script is designed to project the phase mask onto a second monitor (acting as the SLM).

## Example Output

- **Phase Mask**: The generated phase mask is displayed in grayscale. The varying intensities represent the phase shifts needed to shape the beam into the target intensity pattern.
- **Intensity Pattern**: The resulting intensity pattern after applying the phase mask is shown, demonstrating the ability to shape light.

The SLM projection part of this code is based on the `ots1.0.1` toolkit written by Volpe et al.

---

This repository provides a solid foundation for working with optical beam shaping using phase retrieval algorithms and spatial light modulators. Feel free to modify the code for your specific use case, and please cite this repository if you use it in your research or projects.



#### Reference: 
1. https://en.wikipedia.org/wiki/Gerchberg%E2%80%93Saxton_algorithm
2. _Optical Tweezer: Principles and Application_ by Philip H. Jones, O. M. Marago and G. Volpe
   
