# Higher Order Ambisonics Library
#### A compact library for encoding, manipulation and decoding of spatial sound using Higher Order Ambisonics.

---
>
>    Archontis Politis, 2015  
>
>    Department of Signal Processing and Acoustics, Aalto University, Finland  
>
>    archontis.politis@aalto.fi
>
---

This Matlab/Octave library was developed during my doctoral research in the [Communication Acoustics Research Group] (http://spa.aalto.fi/en/research/research_groups/communication_acoustics/), Aalto University, Finland. If you would like to reference the code, you can refer to my dissertation published [here](https://aaltodoc.aalto.fi/handle/123456789/22499):

    Archontis Politis, Microphone array processing for parametric spatial audio techniques, 2016
    Doctoral Dissertation, Department of Signal Processing and Acoustics, Aalto University, Finland
    
## Description

This is a compact Matlab/Octave library implementing most common operations 
associated with higher-order ambisonics (HOA), which refer to a set of
spatial audio techniques for capturing, manipulating and reproducing
sound scenes, based on a spherical Fourier expansion of the sound field.

Most of the functionality of the library is demonstrated in the included 
script TEST_AMBI_SCRIPT.m. The same documentation can be also found in 

[http://research.spa.aalto.fi/projects/ambi-lib/ambi.html]

The included functions implement HOA encoding of directional sounds,
decoding using various decoding approaches, and rotation of HOA sound
scenes. All operations are defined in terms of orthonormalized real 
Spherical Harmonics (N3D in ambisonic slang) and channel indexing
according to $q = n^2+n+m+1$, where n is the order and m is the degree (ACN
in ambisonic slang). However, functions are included to convert to and 
from N3D/ACN to some other established conventions (namely
semi-normalized SHs (SN3D) and an alternative channel indexing, termed
SID).

Ambisonic decoding can be approached from various sides, more physically
inspired or more perceptually inspired. Five approaches are implemented

* 1) Sampling or projection decoding (transpose)
* 2) Mode-matching decoding (pseudo-inverse)
* 3) Energy-preserving decoding [ref.1]
* 4) All-round ambisonic decoding [ref.2]
* 5) Constant-angular spread decoding [ref.3]

Apart from the two first traditional approaches, the three last are more
recent and more perceptually motivated. They are also more
flexible and robust, in terms of loudspeaker layouts.

Additionally, a function evaluating and visualizing the popular ambisonic
performance measures, velocity and energy vectors, along with overall
energy and amplitude preservation, is included.

Max-rE weighting for the decoder [ref.4 & ref.2] can be optionally enabled.

ALLRAD and CSAD decoders require computation of amplitude and energy
panning gains, and large spherical uniform sampling schemes (t-Designs).
Both of these can be found firs in the Matlab/Octave VBAP library in 

[https://github.com/polarch/Vector-Base-Amplitude-Panning]

and the general Spherical harmonic transform library by the author found in

[https://github.com/polarch/Spherical-Harmonic-Transform],

These two libraries should be added to the Matlab path before executing 
this script.

Rotation, apart from the case of simple B-format, also depends on the 
larger spherical harmonic transform library, which contains many 
other operations that may be of interest to ambisonics, like directional 
smoothing (spherical convolution) and directional weighting/shaping 
(spherical multiplication).

The library contains the following main functions:
  
* ambiDecoder:    Compute a HOA decoding matrix for a specified order and
                  a specified method, with or without max-rE weighting
* analyzeDecoder: Analyze amplitude, energy, velocity and energy vector
                  magnitudes and directional errors and spread, for an
                  ambisonic decoder, or from panning gains
* getRSH:         returns values of real orthonormal spherical harmonics
                  vectors of directions
* encodeHOA_N3D:  encode a number of source signals from various directions
                  to arbitrary-order HOA signals (N3D, ACN)
* encodeBformat:  encode a number of source signals from various directions
                  to traditional B-format signals
* decodeHOA_N3D:  decode HOA signals to a loudspeaker setup, using a
                  certain decoding matrix (frequency-dependent decoding
                  possible, see below)
* decodeBformat:  decode B-format signals to a loudspeaker setup, using a
                  certain decoding matrix
* rotateHOA_N3D:  rotate sound scenes encoded or recorded in HOA signals,
                  using a yaw-pitch-roll convention
* rotateBformat:  rotate sound scenes encoded or recorded to B-format signals,
                  using a yaw-pitch-roll convention
* allrad:         implements the all-round ambisonic decoder
* csad:           implements the constant-angular spread decoder
* getLayoutAmbisonicOrder:    Compute the equivalent ambisonic order for
                              regular and irregular speaker arrangements
* getMaxREweights:    Compute the max-rE weights for a certain order, for 
                      decoding (or encoding) weighting
* getTheoreticalEVmag:    Compute the theoretical energy vector magnitude
                          of a plane-wave encoded to a certain order, 
                          with max-rE weighting
* plotSphericalGrid:  Plots angular quantities in a 2D azimuth-elevation
                      grid, with loudspeaker positions superimposed
* convert_ACN_SID:    Convert HOA signals with the ACN indexing to HOA 
                      signals with the SID indexing, and back
* convert_N3D_SN3D:   Convert HOA signals with the N3D normalization to HOA 
                      signals with the SN3D normalization, and back
* convert_N3D_Bformat: Convert 1st-order signals with the N3D and ACN 
                       conventions to traditional B-format, and back


For any questions, comments, corrections, or general feedback, please
contact archontis.politis@aalto.fi


## References

    [1] Zotter, F., Pomberger, H., Noisternig, M. (2012). 
    Energy-Preserving Ambisonic Decoding. Acta Acustica United with Acustica, 98(1), 37:47.

    [2] Zotter, F., Frank, M. (2012). 
    All-Round Ambisonic Panning and Decoding. Journal of the Audio Engineering Society, 60(10), 807:820.

    [3] Epain, N., Jin, C. T., Zotter, F. (2014). 
    Ambisonic Decoding With Constant Angular Spread. Acta Acustica United with Acustica, 100, 928:936.

    [4] Daniel, J. (2001). 
    Representation de champs acoustiques, application ? la transmission et ? la reproduction de scenes sonores complexes dans un contexte multim?dia. Doctoral Thesis. Universite Paris 6.

    [5] Makita, Y. (1962). 
    On the directional localization of sound in the stereophonic sound field. EBU Review, 73, 1536:1539.

    [6] Gerzon, M. A. (1992). 
    General Metatheory of Auditory Localization. In 92nd AES Convention (Preprint 3306). Vienna, Austria.

    [7] Merimaa, J. (2007). 
    Energetic sound field analysis of stereo and multichannel loudspeaker reproduction. In 123rd AES Convention. New York, NY.

    [8] Matthias, F. (2013). 
    Phantom Sources using Multiple Loudspeakers in the Horizontal Plane. Doctoral thesis, Institute of Electronic Music and Acoustics, University of Music and Performing Arts, Graz

    [9] Laitinen, M., Vilkamo, J., Jussila, K., Politis, A., & Pulkki, V. (2014). 
    Gain normalization in amplitude panning as a function of frequency and room reverberance. In 55th International Conference of AES. Helsinki, Finland.

    [10] Frank, M. (2013). 
    Source Width of Frontal Phantom Sources: Perception, Measurement, and Modeling. Archives of Acoustics. 38(3), 311?319

    [11] Daniel, J. (2003). 
    Spatial Sound Encoding Including Near Field Effect : Introducing Distance Coding Filters and a Viable , New Ambisonic Format. In 23rd International Conference of AES. Copenhagen, Denmark.


