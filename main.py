import math
import numpy as np
import matplotlib.pyplot as plt

### GEAR SET 1 ###

# INPUTS

power = 20  # horsepower
lifetime = 12000  # gear lifetime [hours]
inputSpeed = 1750  # rpm
outputSpeed = 85
minNumTeeth = 17  # From table xxxx
Ko = 2.00  # Adjustment factor from 9-1 -> Moderate Shock Input, Moderate Shock Output
serviceFactor = 1.00

# SOLVER

targetVelocityRatio = inputSpeed / outputSpeed
print("Target Gear Ratio: " + str(targetVelocityRatio))

print("GEAR SET 1 PARAMS \n\n\n")

targetVR = math.sqrt(targetVelocityRatio)  # Velocity Ratio for both sets of gears

pinionTeeth = minNumTeeth
gearTeeth = round(pinionTeeth * targetVR)
print("Gear Teeth: " + str(gearTeeth))

VR = gearTeeth / pinionTeeth
print("Velocity Ratio: " + str(VR))

adjustedPower = Ko * power
pitchDiameter = 7.00  # Pd: From 9-11, adjusted power corresponds to a pitch diameter of 7 in

gearDiameter = gearTeeth / pitchDiameter
pinionDiameter = pinionTeeth / pitchDiameter

adjustmentDiameter = 1 / pitchDiameter


print("Adjustment Diameter: " + str(adjustmentDiameter))
print("Gear Diameter: " + str(gearDiameter))
print("Pinion Diameter: " + str(pinionDiameter))

geartrainHeight = (3 / 2) * gearDiameter + pinionDiameter / 2 + 3 * adjustmentDiameter  # calculated Y (in)

print("Geartrain Height: " + str(geartrainHeight))

C = (pinionDiameter + gearDiameter) / (2 * pitchDiameter)  # center distance
tangentialVelocity = (math.pi * pitchDiameter * inputSpeed) / 12
transmittedLoad = 33000 * power / tangentialVelocity  # Wt

print("C: " + str(C))
print("Tangential velocity: " + str(tangentialVelocity))
print("wt: " + str(transmittedLoad))

nominalFaceWidth = 12 / pitchDiameter

print("Nominal Face Width: " + str(nominalFaceWidth))

# Material is steel
elasticCoefficient = 2300  # Cp
agmaQualityNumber = 6  # A6 - automotive transmission equiv. from 9-5

dynamicFactorsAt3207VT = {}  # Mapping between AGMA Quality Numbers and Dynamic Factors at Vt of 3207
dynamicFactorsAt3207VT[6] = 1.13
dynamicFactorsAt3207VT[7] = 1.23
dynamicFactorsAt3207VT[8] = 1.34
dynamicFactorsAt3207VT[9] = 1.46
dynamicFactorsAt3207VT[10] = 1.58
dynamicFactorsAt3207VT[11] = 1.75

dynamicFactor = dynamicFactorsAt3207VT[agmaQualityNumber]  # Kv

print("Dynamic Factor: " + str(dynamicFactor))

KMRatio = nominalFaceWidth / pinionDiameter

print("KM Ratio: " + str(KMRatio))

# Use this formula because 1 < F < 15 (See Figure 9-12)
pinionProportionFactor = nominalFaceWidth / (10 * pinionDiameter) - 0.0375 + 0.0125 * nominalFaceWidth  # from fig 9-12

print("Cpf: " + str(pinionProportionFactor))

# Assume commercial gears
meshAlignmentFactor = 0.127 + 0.0158 * nominalFaceWidth - (1.093 * pow(10, -4) * pow(nominalFaceWidth, 2))  # Cma

Km = 1 + meshAlignmentFactor + pinionProportionFactor

print("Mesh Alignment Factor: " + str(meshAlignmentFactor))

sizeFactor = 1  # Ks since pitch diameter is greater than 5 inches from Table 9-2

rimThicknessFactor = 1.0  # Kb since VR is greater than 1.2, thickness factor is 1.0

reliabilityFactor = 1.00  # Km From table 9-11, at 99% reliability

numPinionCycles = 60 * lifetime * inputSpeed  # Ncp
numGearCycles = 60 * lifetime * (inputSpeed / VR)  # Ncg

print("Number Pinion Cycles: " + str(numPinionCycles))
print("Number Gear Cycles: " + str(numGearCycles))


def bending_strength_stress_cycle_factor(numCycles):
    # This is linear because our load cycles is greater than 10^7 (Fig 9-21)
    return 1.3558 * numCycles ** -0.0178


def pitting_strength_stress_cycle_factor(numCycles):
    # This is linear because our load cycles is greater than 10^7 (Fig 9-22)
    return 1.4488 * numCycles ** -0.023


bendingStressCycleFactorPinion = bending_strength_stress_cycle_factor(numPinionCycles)  # Ynp from Figure 9-21
bendingStressCycleFactorGear = bending_strength_stress_cycle_factor(numGearCycles)  # Yng from Figure 9-21

pittingResistanceStressCycleFactorPinion = pitting_strength_stress_cycle_factor(numPinionCycles)  # Znp from Figure 9-22
pittingResistanceStressCycleFactorGear = pitting_strength_stress_cycle_factor(numGearCycles)  # Zng from Figure 9-22

print("Ynp: " + str(bendingStressCycleFactorPinion))
print("Yng: " + str(bendingStressCycleFactorGear))

print("Znp: " + str(pittingResistanceStressCycleFactorPinion))
print("Zng: " + str(pittingResistanceStressCycleFactorGear))

# From Figure 9-10: Np = 17, Ng = 77
bendingGeometryFactorPinion = 0.295  # Jp from Fig 9-17
bendingGeometryFactorGear = 0.41  # Jg from Fig 9-17

# From Figure 9-17 VR = 4.54 Np = 17

pittingGeometryFactor = 0.105  # I from Figure 9-17, using 20deg pressure angle

allowableBendingStressNumberPinion = (((transmittedLoad * pitchDiameter) / (nominalFaceWidth * bendingGeometryFactorPinion)) * Ko \
                                      * Km * sizeFactor * rimThicknessFactor * dynamicFactor) * (
                                                 serviceFactor * reliabilityFactor / (
                                             bendingStressCycleFactorPinion))  # Satp
allowableBendingStressNumberGear = (((transmittedLoad * pitchDiameter) / (nominalFaceWidth * bendingGeometryFactorGear)) * Ko \
                                    * Km * sizeFactor * rimThicknessFactor * dynamicFactor) * (
                                               serviceFactor * reliabilityFactor / (
                                           bendingStressCycleFactorGear))  # Satg

allowableContactStressNumberPinion = ((elasticCoefficient * reliabilityFactor * serviceFactor) /
                                      pittingResistanceStressCycleFactorPinion) * ((transmittedLoad * Ko * sizeFactor * Km *
                                                                                    dynamicFactor) / (nominalFaceWidth *
                                                                                                      pinionDiameter *
                                                                                                      pittingGeometryFactor)) ** 0.5  # Sacp

allowableContactStressNumberGear = ((elasticCoefficient * reliabilityFactor * serviceFactor) /
                                      pittingResistanceStressCycleFactorGear) * ((transmittedLoad * Ko * sizeFactor * Km *
                                                                                  dynamicFactor) / (nominalFaceWidth *
                                                                                                      gearDiameter *
                                                                                                      pittingGeometryFactor)) ** 0.5  # Sacg

print("Bending Stress Number Pinion: " + str(allowableBendingStressNumberPinion))
print("Bending Stress Number Gear: " + str(allowableBendingStressNumberGear))
print("Contact Stress Number Pinion: " + str(allowableContactStressNumberPinion))
print("Contact Stress Number Gear: " + str(allowableContactStressNumberGear))

limitingStress = max([allowableBendingStressNumberPinion, allowableBendingStressNumberGear, allowableContactStressNumberGear, allowableContactStressNumberPinion])

print("Limiting Factor: " + str(limitingStress))

def calc_brinell_hardness(limitingStress):
    return (limitingStress / 1000 - 29.10) / 0.322

brinellHardness = calc_brinell_hardness(limitingStress) # HB from Fig 9-19 using Grade 1 Steel
print("Brinell Hardness: " + str(brinellHardness))

# Use appendix 3 to select a steel -> Use 1144 Cold-Drawn (HB = 200)


### GEAR SET 2 ###

print("\n\n\n GEAR SET 2 PARAMS \n\n\n ")

targetVR = math.sqrt(targetVelocityRatio)  # Velocity Ratio for both sets of gears

pinionTeeth = minNumTeeth
gearTeeth = round(pinionTeeth * targetVR)
print("Gear Teeth: " + str(gearTeeth))

VR = gearTeeth / pinionTeeth
print("Velocity Ratio: " + str(VR))

inputSpeed = inputSpeed / VR

adjustedPower = Ko * power
pitchDiameter = 7.00  # Pd: From 9-11, adjusted power corresponds to a pitch diameter of 7 in

gearDiameter = gearTeeth / pitchDiameter
pinionDiameter = pinionTeeth / pitchDiameter

adjustmentDiameter = 1 / pitchDiameter

print("Gear Diameter: " + str(gearDiameter))
print("Pinion Diameter: " + str(pinionDiameter))

geartrainHeight = (3 / 2) * gearDiameter + pinionDiameter / 2 + 3 * adjustmentDiameter  # calculated Y (in)

print("Geartrain Height: " + str(geartrainHeight))

C = (pinionDiameter + gearDiameter) / (2 * pitchDiameter)  # center distance
tangentialVelocity = (math.pi * pitchDiameter * inputSpeed) / 12
transmittedLoad = 33000 * power / tangentialVelocity  # Wt

print("C: " + str(C))
print("Tangential velocity: " + str(tangentialVelocity))
print("wt: " + str(transmittedLoad))

nominalFaceWidth = 12 / pitchDiameter

print("Nominal Face Width: " + str(nominalFaceWidth))

# Material is steel
elasticCoefficient = 2300  # Cp
agmaQualityNumber = 6  # A6 - automotive transmission equiv. from 9-5

dynamicFactorsAt3207VT = {}  # Mapping between AGMA Quality Numbers and Dynamic Factors at Vt of 708
dynamicFactorsAt3207VT[6] = 1.07
dynamicFactorsAt3207VT[7] = 1.12
dynamicFactorsAt3207VT[8] = 1.18
dynamicFactorsAt3207VT[9] = 1.24
dynamicFactorsAt3207VT[10] = 1.31
dynamicFactorsAt3207VT[11] = 1.39

dynamicFactor = dynamicFactorsAt3207VT[agmaQualityNumber]  # Kv

print("Dynamic Factor: " + str(dynamicFactor))

KMRatio = nominalFaceWidth / pinionDiameter

print("KM Ratio: " + str(KMRatio))

# Use this formula because 1 < F < 15 (See Figure 9-12)
pinionProportionFactor = nominalFaceWidth / (10 * pinionDiameter) - 0.0375 + 0.0125 * nominalFaceWidth  # from fig 9-12

print("Cpf: " + str(pinionProportionFactor))

# Assume commercial gears
meshAlignmentFactor = 0.127 + 0.0158 * nominalFaceWidth - (1.093 * pow(10, -4) * pow(nominalFaceWidth, 2))  # Cma

Km = 1 + meshAlignmentFactor + pinionProportionFactor

print("Mesh Alignment Factor: " + str(meshAlignmentFactor))

sizeFactor = 1  # Ks since pitch diameter is greater than 5 inches from Table 9-2

rimThicknessFactor = 1.0  # Kb since VR is greater than 1.2, thickness factor is 1.0

reliabilityFactor = 1.00  # Km From table 9-11, at 99% reliability

numPinionCycles = 60 * lifetime * inputSpeed  # Ncp
numGearCycles = 60 * lifetime * (inputSpeed / VR)  # Ncg

print("Number Pinion Cycles: " + str(numPinionCycles))
print("Number Gear Cycles: " + str(numGearCycles))

bendingStressCycleFactorPinion = bending_strength_stress_cycle_factor(numPinionCycles)  # Ynp from Figure 9-21
bendingStressCycleFactorGear = bending_strength_stress_cycle_factor(numGearCycles)  # Yng from Figure 9-21

pittingResistanceStressCycleFactorPinion = pitting_strength_stress_cycle_factor(numPinionCycles)  # Znp from Figure 9-22
pittingResistanceStressCycleFactorGear = pitting_strength_stress_cycle_factor(numGearCycles)  # Zng from Figure 9-22

print("Ynp: " + str(bendingStressCycleFactorPinion))
print("Yng: " + str(bendingStressCycleFactorGear))

print("Znp: " + str(pittingResistanceStressCycleFactorPinion))
print("Zng: " + str(pittingResistanceStressCycleFactorGear))

# From Figure 9-10: Np = 17, Ng = 77
bendingGeometryFactorPinion = 0.295  # Jp from Fig 9-17
bendingGeometryFactorGear = 0.41  # Jg from Fig 9-17

# From Figure 9-17 VR = 4.54 Np = 17

pittingGeometryFactor = 0.105  # I from Figure 9-17, using 20deg pressure angle

allowableBendingStressNumberPinion = (((transmittedLoad * pitchDiameter) / (nominalFaceWidth * bendingGeometryFactorPinion)) * Ko \
                                      * Km * sizeFactor * rimThicknessFactor * dynamicFactor) * (
                                                 serviceFactor * reliabilityFactor / (
                                             bendingStressCycleFactorPinion))  # Satp
allowableBendingStressNumberGear = (((transmittedLoad * pitchDiameter) / (nominalFaceWidth * bendingGeometryFactorGear)) * Ko \
                                    * Km * sizeFactor * rimThicknessFactor * dynamicFactor) * (
                                               serviceFactor * reliabilityFactor / (
                                           bendingStressCycleFactorGear))  # Satg

allowableContactStressNumberPinion = ((elasticCoefficient * reliabilityFactor * serviceFactor) /
                                      pittingResistanceStressCycleFactorPinion) * ((transmittedLoad * Ko * sizeFactor * Km *
                                                                                    dynamicFactor) / (nominalFaceWidth *
                                                                                                      pinionDiameter *
                                                                                                      pittingGeometryFactor)) ** 0.5  # Sacp

allowableContactStressNumberGear = ((elasticCoefficient * reliabilityFactor * serviceFactor) /
                                      pittingResistanceStressCycleFactorGear) * ((transmittedLoad * Ko * sizeFactor * Km *
                                                                                  dynamicFactor) / (nominalFaceWidth *
                                                                                                      gearDiameter *
                                                                                                      pittingGeometryFactor)) ** 0.5  # Sacg

print("Bending Stress Number Pinion: " + str(allowableBendingStressNumberPinion))
print("Bending Stress Number Gear: " + str(allowableBendingStressNumberGear))
print("Contact Stress Number Pinion: " + str(allowableContactStressNumberPinion))
print("Contact Stress Number Gear: " + str(allowableContactStressNumberGear))

limitingStress = max([allowableBendingStressNumberPinion, allowableBendingStressNumberGear, allowableContactStressNumberGear, allowableContactStressNumberPinion])

print("Limiting Factor: " + str(limitingStress))

brinellHardness = calc_brinell_hardness(limitingStress) # HB from Fig 9-19 using Grade 1 Steel
print("Brinell Hardness: " + str(brinellHardness))

print("Output Speed: " + str(inputSpeed / VR))
# Use appendix 3 to select a steel -> Use SAE 4150 OQT 700 (HB = 495)

# Shaft 1

print("\n\n\nSHAFT 1\n\n\n")

powerIn = 20 * 63000 # watts
lengthSec1 = 0.5 # inches
lengthSec2 = 1
lengthSec3ToMiddle = 1.107
lengthSec3MiddleToEnd = 0.857
lengthSec4 = 1
lengthSec5 = 0.25
lengthSec6 = 0.25

lengthAToCenter = lengthSec1 + lengthSec2 + lengthSec3ToMiddle
lengthCenterToC = lengthSec3MiddleToEnd + lengthSec4 + lengthSec5 + lengthSec6;
totalLength = lengthAToCenter + lengthCenterToC

def plot_shear_diagram(tangentialForce, lengthAToCenter, lengthCenterToC, figureNumber, plotTitle):
    rc = tangentialForce * lengthAToCenter / (lengthAToCenter + lengthCenterToC)
    ra = tangentialForce - rc
    x = np.linspace(0, lengthAToCenter + lengthCenterToC, discreteSubdivisions)
    y = np.linspace(0, 0, discreteSubdivisions)
    for i in range(0, int(discreteSubdivisions * lengthAToCenter / (lengthAToCenter + lengthCenterToC))) :
        y[i] = ra

    for i in range(int(discreteSubdivisions * lengthAToCenter / (lengthAToCenter + lengthCenterToC)), discreteSubdivisions):
        y[i] = -1 * rc

    y[499] = 0
    y[0] = 0
    plt.figure(figureNumber)
    plt.title(plotTitle)
    plt.ylabel("Shear Force [lbs]")
    plt.xlabel("x [in]")
    plt.plot(x, y)

    return [x, y]

def plot_moment_diagram(distance, shearForce, figureNumber, plotTitle):
    moment = np.linspace(0, 0, discreteSubdivisions)
    for i in range(0, discreteSubdivisions):
        moment[i] = np.trapz(shearForce[0 : i], distance[0 : i])
    plt.figure(figureNumber)
    plt.title(plotTitle)
    plt.ylabel("Moment [lbs-in]")
    plt.xlabel("x [in]")
    plt.plot(distance, moment)

    # Adjust by Kt
    for i in range(int(discreteSubdivisions * (lengthSec1 + lengthSec2) / (lengthAToCenter + lengthCenterToC)),
                   int(discreteSubdivisions * (lengthAToCenter + lengthSec3MiddleToEnd) / (lengthAToCenter + lengthCenterToC))):
        moment[i] = moment[i] * Kt
    plt.figure(figureNumber + 1)
    plt.title(plotTitle + " - Kt adjusted")
    plt.ylabel("Moment [lbs-in]")
    plt.xlabel("x [in]")
    plt.plot(distance, moment)
    return [distance, moment]

def cs_from_d(D): # From Table 5-4
    if D <= 0.3:
        return 1
    elif D <= 2 :
        return (D/0.3)**-0.11
    else:
        return 0.859 - 0.02125 * D


def get_section_diameter(torque, moment):
    return ((32 * N / math.pi) * ((Kt * moment / SnPrime) ** 2 + 0.75 * (torque / Sy)**2) ** 1/2) ** (1/3)

def get_section_diameters(tangentialMoment, radialMoment, torque):
    section1MaxMoment = 0
    for i in range(0, int(discreteSubdivisions * lengthSec1 / totalLength)):
        if math.sqrt(tangentialMoment[i] ** 2 + radialMoment[i] ** 2) > section1MaxMoment:
            section1MaxMoment = math.sqrt(tangentialMoment[i] ** 2 + radialMoment[i] ** 2)
    section1Diameter = get_section_diameter(torqueInB, section1MaxMoment)

    section2MaxMoment = 0
    for i in range(int(discreteSubdivisions * lengthSec1 / totalLength),
                   int(discreteSubdivisions * (lengthSec1 + lengthSec2) / totalLength)):
        if math.sqrt(tangentialMoment[i] ** 2 + radialMoment[i] ** 2) > section2MaxMoment:
            section2MaxMoment = math.sqrt(tangentialMoment[i] ** 2 + radialMoment[i] ** 2)
    section2Diameter = get_section_diameter(torqueInB, section2MaxMoment)

    section3MaxMoment = 0

    for i in range(int(discreteSubdivisions * (lengthSec1 + lengthSec2) / totalLength),
                   int(discreteSubdivisions * (lengthSec1 + lengthSec2 + lengthSec3MiddleToEnd + lengthSec3ToMiddle) / totalLength)):
        if math.sqrt(tangentialMoment[i] ** 2 + radialMoment[i] ** 2) > section3MaxMoment:
            section3MaxMoment = math.sqrt(tangentialMoment[i] ** 2 + radialMoment[i] ** 2)
    section3Diameter = get_section_diameter(torqueInB, section3MaxMoment)

    section4MaxMoment = 0
    for i in range(int(discreteSubdivisions * (lengthSec1 + lengthSec2 + lengthSec3MiddleToEnd + lengthSec3ToMiddle) / totalLength),
                   int(discreteSubdivisions * (lengthSec1 + lengthSec2 + lengthSec3MiddleToEnd + lengthSec3ToMiddle + lengthSec4) / totalLength)):
        if math.sqrt(tangentialMoment[i] ** 2 + radialMoment[i] ** 2) > section4MaxMoment:
            section4MaxMoment = math.sqrt(tangentialMoment[i] ** 2 + radialMoment[i] ** 2)
    section4Diameter = get_section_diameter(0, section4MaxMoment)

    section5MaxMoment = 0
    for i in range(int(discreteSubdivisions * (lengthSec1 + lengthSec2 + lengthSec3MiddleToEnd + lengthSec3ToMiddle + lengthSec4) / totalLength),
                   int(discreteSubdivisions * (lengthSec1 + lengthSec2 + lengthSec3MiddleToEnd + lengthSec3ToMiddle + lengthSec4 + lengthSec5) / totalLength)):
        if math.sqrt(tangentialMoment[i] ** 2 + radialMoment[i] ** 2) > section5MaxMoment:
            section5MaxMoment = math.sqrt(tangentialMoment[i] ** 2 + radialMoment[i] ** 2)
    section5Diameter = get_section_diameter(0, section5MaxMoment)


    section6MaxMoment = 0
    for i in range(int(discreteSubdivisions * (lengthSec1 + lengthSec2 + lengthSec3MiddleToEnd + lengthSec3ToMiddle + lengthSec4 + lengthSec5) / totalLength),
                   int(discreteSubdivisions)):
        if math.sqrt(tangentialMoment[i] ** 2 + radialMoment[i] ** 2) > section6MaxMoment:
            section6MaxMoment = math.sqrt(tangentialMoment[i] ** 2 + radialMoment[i] ** 2)
    section6Diameter = get_section_diameter(0, section6MaxMoment)

    return [section1Diameter, section2Diameter, section3Diameter, section4Diameter, section5Diameter, section6Diameter]

inputSpeed = 1750
Kt = 3
Cr = 0.81 # 99% reliability from Table 5-3
D = 2
Sn = 80000 # SAE 1040 Cold-Drawn (Appendix 3) [ksi]
Sy = 71000 # SAE 1040 Cold-Drawn (Appendix 3) [ksi]
discreteSubdivisions = 500

Cs = cs_from_d(D)
SnPrime = Cs * Cr * Sn
N = 3 # Safety Factor

angularVelocityIn = inputSpeed * 2 * math.pi / 60 # rpm
torqueInB = powerIn / angularVelocityIn
print("Torque In: " + str(torqueInB))

tangentialForceGear = torqueInB / (pinionDiameter / 2)
radialForceGear = tangentialForceGear * math.tan(math.pi / 9)

print("Tangential Force Gear: " + str(tangentialForceGear))
print("Radial Force Gear: " + str(radialForceGear))

def father_of_get_section_diameters(tangentialForceGear, radialForceGear, startingFigureNumber, torque):
    [distance, shearForceTangential] = plot_shear_diagram(tangentialForceGear, lengthAToCenter, lengthCenterToC, startingFigureNumber, "Tangential Shear Force")
    [distance, shearForceRadial] = plot_shear_diagram(radialForceGear, lengthAToCenter, lengthCenterToC, startingFigureNumber + 1, "Radial Shear Force")

    [distance, tangentialMoment] = plot_moment_diagram(distance, shearForceTangential, startingFigureNumber + 2, "Moment Diagram Tangential")
    [distance, radialMoment] = plot_moment_diagram(distance, shearForceRadial, startingFigureNumber + 4, "Moment Diagram Radial")

    return get_section_diameters(tangentialMoment, radialMoment, torque)

section_diameters = father_of_get_section_diameters(tangentialForceGear, radialForceGear, 1, torqueInB)
print("Input Section Diameters: " + str(section_diameters))

angularVelocityIn = (inputSpeed / VR) * 2 * math.pi / 60
torqueInB = powerIn / angularVelocityIn
print("Torque In: " + str(torqueInB))

tangentialForceGear = torqueInB / (gearDiameter / 2)
radialForceGear = tangentialForceGear * math.tan(math.pi / 9)
print("Tangential Force Gear: " + str(tangentialForceGear))
print("Radial Force Gear: " + str(radialForceGear))

section_diameters = father_of_get_section_diameters(tangentialForceGear, radialForceGear, 7, torqueInB)
print("Output Section Diameters: " + str(section_diameters))

plt.show()
