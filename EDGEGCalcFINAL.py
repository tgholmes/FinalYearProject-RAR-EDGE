# This code calculates the baryonic and observed gravitational accelerations on the reduced HiRes files, that can be then exported for further analysis


import numpy as np

def CalculateGbaryonic(galaxy):
    starsData = np.loadtxt(f'{galaxy}/{galaxy}_star_mass_HIRES.txt')
    gasData = np.loadtxt(f'{galaxy}/{galaxy}_gas_mass_HIRES.txt')
    
    starsRadii = starsData[:, 0]
    starsMass = starsData[:, 1]
    gasRadii = gasData[:, 0]
    gasMass = gasData[:, 1]
    
    G = 6.67e-11
    rSteps = np.logspace(np.log10(100), np.log10(5000), 25)
    
    CombinedR = np.concatenate((starsRadii, gasRadii))
    CombinedMass = np.concatenate((starsMass, gasMass))
    
    SortedIndices = np.argsort(CombinedR)
    SortedR = CombinedR[SortedIndices]
    SortedMass = CombinedMass[SortedIndices]
    
    cumulativeMassBaryonic = np.zeros_like(rSteps)
    gbar = np.zeros_like(rSteps)
    for i in range(len(rSteps)):
        particles_within_step = SortedMass[SortedR <= rSteps[i]]
        cumulativeMassBaryonic[i] = np.sum(particles_within_step)
        gbar[i] = (G * cumulativeMassBaryonic[i] * 1.989e30) / (rSteps[i] * 3.086e16)**2

    return gbar, rSteps

#---------------------------------------------------------------------------------

# Similar function to calculate g observed
def CalculateGobserved(galaxy):
    starsData = np.loadtxt(f'{galaxy}/{galaxy}_star_mass_HIRES.txt')
    gasData = np.loadtxt(f'{galaxy}/{galaxy}_gas_mass_HIRES.txt')
    dmData = np.loadtxt(f'{galaxy}/{galaxy}_dm_mass_HIRES.txt')
    
    starsRadii = starsData[:, 0]
    starsMass = starsData[:, 1]
    gasRadii = gasData[:, 0]
    gasMass = gasData[:, 1]
    dmRadii = dmData[:, 0]
    dmMass = dmData[:, 1]
    
    G = 6.67e-11
    rSteps = np.logspace(np.log10(100), np.log10(5000), 25)

    CombinedR = np.concatenate((starsRadii, gasRadii, dmRadii))
    CombinedMass = np.concatenate((starsMass, gasMass, dmMass))
    SortedIndices = np.argsort(CombinedR)
    SortedR = CombinedR[SortedIndices]
    SortedMass = CombinedMass[SortedIndices]
    
    cumulativeMass = np.zeros_like(rSteps)
    gobs = np.zeros_like(rSteps)
    for i in range(len(rSteps)):
        particles_within_step = SortedMass[SortedR <= rSteps[i]]
        cumulativeMass[i] = np.sum(particles_within_step) #Use indexing instead
        gobs[i] = (G * cumulativeMass[i] * 1.989e30) / (rSteps[i] * 3.086e16)**2
    
    return gobs, rSteps

galaxy = '1459'

gbar, r_steps_gbar = CalculateGbaryonic(galaxy)
gobs, r_steps_gobs = CalculateGobserved(galaxy)

# Write data to file
output_file = f'{galaxy}/{galaxy}_g_data_HIRES_TRIAL.txt'
np.savetxt(output_file, np.column_stack((r_steps_gbar, gbar, gobs)), fmt='%.18e')

