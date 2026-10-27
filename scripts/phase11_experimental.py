#!/usr/bin/env python3
"""
Phase 11: Experimental Apparatus Design
========================================
Instrumentation design, safety protocols, and data acquisition for
backpack-scale spacetime manipulation experiment.

No human subjects - instrument-only validation.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional
import json
import math


@dataclass
class SensorSpec:
    """Sensor specification for experimental apparatus."""
    name: str
    quantity: str
    range_min: float
    range_max: float
    accuracy: float
    sample_rate_hz: float
    interface: str
    power_w: float
    mass_g: float
    cost_usd: float
    notes: str = ""


@dataclass
class SafetyInterlock:
    """Safety interlock specification."""
    name: str
    trigger_condition: str
    response_time_ms: float
    action: str
    redundancy_level: int  # 1=single, 2=dual, 3=triple
    notes: str = ""


@dataclass
class DAQChannel:
    """Data acquisition channel specification."""
    name: str
    sensor: str
    channel_type: str  # analog_in, digital_in, analog_out, digital_out
    sample_rate_hz: float
    resolution_bits: int
    voltage_range_v: tuple
    notes: str = ""


@dataclass
class ExperimentalApparatus:
    """Complete experimental apparatus specification."""
    # Coil system (from Phase 10)
    coil_spec: Dict
    
    # Power system
    hv_supply_spec: Dict
    battery_spec: Dict
    pfc_capacitor_f: float
    resonant_capacitor_f: float
    
    # Sensors
    magnetic_sensors: List[SensorSpec]
    thermal_sensors: List[SensorSpec]
    electrical_sensors: List[SensorSpec]
    mechanical_sensors: List[SensorSpec]
    environmental_sensors: List[SensorSpec]
    
    # Safety
    interlocks: List[SafetyInterlock]
    emergency_stop: SafetyInterlock
    
    # Data acquisition
    daq_channels: List[DAQChannel]
    daq_sample_rate_hz: float
    daq_buffer_size: int
    storage_gb_per_hour: float
    
    # Mechanical
    enclosure_spec: Dict
    cooling_spec: Dict
    total_mass_kg: float
    total_volume_l: float
    
    # Test procedures
    test_procedures: List[Dict]
    
    # Compliance
    safety_standards: List[str]
    emc_standards: List[str]


def design_magnetic_sensors() -> List[SensorSpec]:
    """Design magnetic field sensor suite."""
    return [
        SensorSpec(
            name="Fluxgate Magnetometer (3-axis)",
            quantity="B-field vector",
            range_min=-100e-6,  # -100 uT
            range_max=100e-6,   # +100 uT
            accuracy=10e-9,     # 10 nT
            sample_rate_hz=10000,
            interface="SPI",
            power_w=0.15,
            mass_g=25,
            cost_usd=180,
            notes="Primary B-field measurement at coil center and gradient positions"
        ),
        SensorSpec(
            name="Hall Effect Sensor Array (8x)",
            quantity="B-field mapping",
            range_min=-50e-3,   # -50 mT
            range_max=50e-3,    # +50 mT
            accuracy=100e-6,    # 100 uT
            sample_rate_hz=50000,
            interface="I2C",
            power_w=0.05,
            mass_g=5,
            cost_usd=12,
            notes="Spatial mapping around coil; 8 sensors in 3D grid"
        ),
        SensorSpec(
            name="Search Coil (dB/dt)",
            quantity="dB/dt",
            range_min=-10,
            range_max=10,
            accuracy=0.01,
            sample_rate_hz=200000,
            interface="BNC/ADC",
            power_w=0.0,
            mass_g=50,
            cost_usd=85,
            notes="High-bandwidth dB/dt for pulse characterization"
        ),
        SensorSpec(
            name="NV Center Magnetometer (optional)",
            quantity="B-field vector",
            range_min=-1e-3,
            range_max=1e-3,
            accuracy=1e-12,     # 1 pT/rtHz
            sample_rate_hz=1000,
            interface="USB",
            power_w=2.0,
            mass_g=200,
            cost_usd=15000,
            notes="Quantum sensor for ultra-low field validation; optional upgrade"
        ),
    ]


def design_thermal_sensors() -> List[SensorSpec]:
    """Design thermal sensor suite."""
    return [
        SensorSpec(
            name="RTD Pt100 (4-wire) x12",
            quantity="Temperature",
            range_min=-50,
            range_max=200,
            accuracy=0.1,
            sample_rate_hz=100,
            interface="SPI/ADC",
            power_w=0.01,
            mass_g=10,
            cost_usd=25,
            notes="Coil winding (4), core (2), power electronics (4), ambient (2)"
        ),
        SensorSpec(
            name="Thermocouple Type K x8",
            quantity="Temperature",
            range_min=-200,
            range_max=1200,
            accuracy=1.5,
            sample_rate_hz=1000,
            interface="ADC",
            power_w=0.0,
            mass_g=5,
            cost_usd=8,
            notes="High-temp backup on coil surface and power components"
        ),
        SensorSpec(
            name="IR Thermal Camera (32x32)",
            quantity="Surface temp map",
            range_min=-20,
            range_max=300,
            accuracy=1.0,
            sample_rate_hz=9,
            interface="I2C",
            power_w=0.5,
            mass_g=45,
            cost_usd=120,
            notes="Non-contact surface temperature mapping"
        ),
    ]


def design_electrical_sensors() -> List[SensorSpec]:
    """Design electrical measurement sensor suite."""
    return [
        SensorSpec(
            name="High-Voltage Differential Probe",
            quantity="Voltage",
            range_min=-1000,
            range_max=1000,
            accuracy=0.5,
            sample_rate_hz=500000,
            interface="BNC/ADC",
            power_w=0.5,
            mass_g=150,
            cost_usd=450,
            notes="Coil voltage waveform capture"
        ),
        SensorSpec(
            name="Current Transducer (Hall, 200A)",
            quantity="Current",
            range_min=-200,
            range_max=200,
            accuracy=0.5,
            sample_rate_hz=200000,
            interface="BNC/ADC",
            power_w=0.3,
            mass_g=80,
            cost_usd=180,
            notes="Coil current waveform; galvanic isolation"
        ),
        SensorSpec(
            name="Rogowski Coil (500A)",
            quantity="Current (AC)",
            range_min=-500,
            range_max=500,
            accuracy=1.0,
            sample_rate_hz=1000000,
            interface="BNC/ADC",
            power_w=0.1,
            mass_g=60,
            cost_usd=220,
            notes="High-bandwidth AC current for harmonic analysis"
        ),
        SensorSpec(
            name="Power Analyzer (3-phase)",
            quantity="P, Q, S, PF",
            range_min=0,
            range_max=5000,
            accuracy=0.2,
            sample_rate_hz=10000,
            interface="Ethernet",
            power_w=5.0,
            mass_g=800,
            cost_usd=3500,
            notes="Real-time power quality monitoring"
        ),
        SensorSpec(
            name="Battery Monitor (Coulomb Counter)",
            quantity="SoC, Current, Voltage",
            range_min=0,
            range_max=100,
            accuracy=1.0,
            sample_rate_hz=10,
            interface="I2C",
            power_w=0.02,
            mass_g=15,
            cost_usd=35,
            notes="Battery state of charge and health tracking"
        ),
    ]


def design_mechanical_sensors() -> List[SensorSpec]:
    """Design mechanical/vibration sensor suite."""
    return [
        SensorSpec(
            name="MEMS Accelerometer (3-axis) x4",
            quantity="Acceleration",
            range_min=-16,
            range_max=16,
            accuracy=0.001,
            sample_rate_hz=10000,
            interface="SPI",
            power_w=0.01,
            mass_g=3,
            cost_usd=15,
            notes="Coil vibration (2), enclosure (2); structural monitoring"
        ),
        SensorSpec(
            name="Strain Gauge (full bridge) x6",
            quantity="Strain",
            range_min=-5000,
            range_max=5000,
            accuracy=1.0,
            sample_rate_hz=10000,
            interface="ADC",
            power_w=0.05,
            mass_g=2,
            cost_usd=45,
            notes="Coil former strain, mounting points"
        ),
        SensorSpec(
            name="Laser Displacement Sensor",
            quantity="Displacement",
            range_min=-5,
            range_max=5,
            accuracy=1e-6,
            sample_rate_hz=50000,
            interface="Ethernet",
            power_w=2.0,
            mass_g=300,
            cost_usd=2500,
            notes="Coil radial expansion under Lorentz forces"
        ),
    ]


def design_environmental_sensors() -> List[SensorSpec]:
    """Design environmental sensor suite."""
    return [
        SensorSpec(
            name="Barometric Pressure Sensor",
            quantity="Pressure",
            range_min=300,
            range_max=1100,
            accuracy=0.5,
            sample_rate_hz=10,
            interface="I2C",
            power_w=0.003,
            mass_g=2,
            cost_usd=8,
            notes="Altitude/pressure compensation for magnetic sensors"
        ),
        SensorSpec(
            name="Humidity/Temperature Sensor",
            quantity="RH, Temp",
            range_min=0,
            range_max=100,
            accuracy=2.0,
            sample_rate_hz=1,
            interface="I2C",
            power_w=0.001,
            mass_g=3,
            cost_usd=6,
            notes="Environmental monitoring for sensor drift correction"
        ),
        SensorSpec(
            name="3-Axis Magnetometer (Earth field)",
            quantity="Earth B-field",
            range_min=-50e-6,
            range_max=50e-6,
            accuracy=50e-9,
            sample_rate_hz=100,
            interface="I2C",
            power_w=0.002,
            mass_g=5,
            cost_usd=12,
            notes="Earth field reference for background subtraction"
        ),
    ]


def design_safety_interlocks() -> List[SafetyInterlock]:
    """Design safety interlock system."""
    return [
        SafetyInterlock(
            name="Coil Overcurrent",
            trigger_condition="I_coil > 120 A (120% rated)",
            response_time_ms=0.1,
            action="Immediate HV supply shutdown + crowbar",
            redundancy_level=3,
            notes="Triple-redundant: HW comparator + FW watchdog + SW monitor"
        ),
        SafetyInterlock(
            name="Coil Overtemperature",
            trigger_condition="T_winding > 145 C (Class F limit - 10 C margin)",
            response_time_ms=100,
            action="Reduce duty cycle to 0% + forced air max",
            redundancy_level=2,
            notes="Dual RTD per winding location; voting logic"
        ),
        SafetyInterlock(
            name="Power Supply Overvoltage",
            trigger_condition="V_coil > 60 V (120% rated)",
            response_time_ms=0.05,
            action="HV supply shutdown + crowbar",
            redundancy_level=3,
            notes="Fast HW overvoltage protection"
        ),
        SafetyInterlock(
            name="Battery Undervoltage",
            trigger_condition="V_batt < 3.0 V/cell",
            response_time_ms=1000,
            action="Graceful shutdown sequence",
            redundancy_level=2,
            notes="Prevents deep discharge damage"
        ),
        SafetyInterlock(
            name="Battery Overtemperature",
            trigger_condition="T_batt > 55 C",
            response_time_ms=500,
            action="Reduce power + increase cooling",
            redundancy_level=2,
            notes="Li-ion thermal runaway prevention"
        ),
        SafetyInterlock(
            name="Coolant Flow Loss (if liquid)",
            trigger_condition="Flow < 0.5 L/min",
            response_time_ms=200,
            action="Immediate power reduction to 10%",
            redundancy_level=2,
            notes="Flow sensor + temperature rise cross-check"
        ),
        SafetyInterlock(
            name="Enclosure Pressure Rise",
            trigger_condition="dP/dt > 1 kPa/s",
            response_time_ms=50,
            action="Emergency vent + power cutoff",
            redundancy_level=2,
            notes="Arc fault / catastrophic failure detection"
        ),
        SafetyInterlock(
            name="EMI/RFI Threshold",
            trigger_condition="Radiated > FCC Part 15 Class A limit",
            response_time_ms=10,
            action="Frequency hop / power reduction",
            redundancy_level=1,
            notes="Real-time spectrum monitor; adaptive mitigation"
        ),
    ]


def design_emergency_stop() -> SafetyInterlock:
    """Design emergency stop system."""
    return SafetyInterlock(
        name="Emergency Stop (E-Stop)",
        trigger_condition="Manual button press OR remote wireless signal OR watchdog timeout",
        response_time_ms=1,
        action="Hard power cutoff (contactor open) + capacitor discharge",
        redundancy_level=3,
        notes="Hardwired, failsafe, latching; requires manual reset; wireless backup"
    )


def design_daq_channels() -> List[DAQChannel]:
    """Design data acquisition channels."""
    channels = []
    
    # Magnetic sensors
    channels.append(DAQChannel(
        name="Fluxgate_X", sensor="Fluxgate Magnetometer", channel_type="analog_in",
        sample_rate_hz=10000, resolution_bits=24, voltage_range_v=(-10, 10),
        notes="Bx component"
    ))
    channels.append(DAQChannel(
        name="Fluxgate_Y", sensor="Fluxgate Magnetometer", channel_type="analog_in",
        sample_rate_hz=10000, resolution_bits=24, voltage_range_v=(-10, 10),
        notes="By component"
    ))
    channels.append(DAQChannel(
        name="Fluxgate_Z", sensor="Fluxgate Magnetometer", channel_type="analog_in",
        sample_rate_hz=10000, resolution_bits=24, voltage_range_v=(-10, 10),
        notes="Bz component"
    ))
    
    # Hall array (8 sensors x 3 axes = 24 channels, multiplexed)
    for i in range(8):
        for axis in ['X', 'Y', 'Z']:
            channels.append(DAQChannel(
                name=f"Hall_{i}_{axis}", sensor=f"Hall Effect Sensor {i}", channel_type="analog_in",
                sample_rate_hz=50000, resolution_bits=16, voltage_range_v=(-5, 5),
                notes=f"Hall sensor {i} {axis}-axis"
            ))
    
    # Search coil
    channels.append(DAQChannel(
        name="Search_Coil", sensor="Search Coil", channel_type="analog_in",
        sample_rate_hz=200000, resolution_bits=16, voltage_range_v=(-10, 10),
        notes="dB/dt high bandwidth"
    ))
    
    # Electrical
    channels.append(DAQChannel(
        name="Coil_Voltage", sensor="HV Differential Probe", channel_type="analog_in",
        sample_rate_hz=500000, resolution_bits=16, voltage_range_v=(-10, 10),
        notes="Coil voltage waveform"
    ))
    channels.append(DAQChannel(
        name="Coil_Current_Hall", sensor="Current Transducer", channel_type="analog_in",
        sample_rate_hz=200000, resolution_bits=16, voltage_range_v=(-10, 10),
        notes="Coil current (Hall)"
    ))
    channels.append(DAQChannel(
        name="Coil_Current_Rogowski", sensor="Rogowski Coil", channel_type="analog_in",
        sample_rate_hz=1000000, resolution_bits=16, voltage_range_v=(-10, 10),
        notes="Coil current (Rogowski, high BW)"
    ))
    
    # Thermal (RTDs - 12 channels)
    for i in range(12):
        channels.append(DAQChannel(
            name=f"RTD_{i}", sensor="RTD Pt100", channel_type="analog_in",
            sample_rate_hz=100, resolution_bits=24, voltage_range_v=(0, 5),
            notes=f"RTD channel {i}"
        ))
    
    # Thermocouples (8 channels)
    for i in range(8):
        channels.append(DAQChannel(
            name=f"TC_{i}", sensor="Thermocouple Type K", channel_type="analog_in",
            sample_rate_hz=1000, resolution_bits=16, voltage_range_v=(-0.05, 0.05),
            notes=f"Thermocouple {i}"
        ))
    
    # Mechanical
    for i in range(4):
        for axis in ['X', 'Y', 'Z']:
            channels.append(DAQChannel(
                name=f"Accel_{i}_{axis}", sensor=f"MEMS Accelerometer {i}", channel_type="analog_in",
                sample_rate_hz=10000, resolution_bits=16, voltage_range_v=(-5, 5),
                notes=f"Accelerometer {i} {axis}"
            ))
    
    for i in range(6):
        channels.append(DAQChannel(
            name=f"Strain_{i}", sensor=f"Strain Gauge {i}", channel_type="analog_in",
            sample_rate_hz=10000, resolution_bits=24, voltage_range_v=(-0.05, 0.05),
            notes=f"Strain gauge {i}"
        ))
    
    channels.append(DAQChannel(
        name="Laser_Displacement", sensor="Laser Displacement Sensor", channel_type="analog_in",
        sample_rate_hz=50000, resolution_bits=16, voltage_range_v=(-10, 10),
        notes="Coil radial displacement"
    ))
    
    # Environmental
    channels.append(DAQChannel(
        name="Baro_Pressure", sensor="Barometric Pressure", channel_type="analog_in",
        sample_rate_hz=10, resolution_bits=16, voltage_range_v=(0, 5),
        notes="Atmospheric pressure"
    ))
    channels.append(DAQChannel(
        name="Humidity", sensor="Humidity Sensor", channel_type="analog_in",
        sample_rate_hz=1, resolution_bits=12, voltage_range_v=(0, 3.3),
        notes="Relative humidity"
    ))
    channels.append(DAQChannel(
        name="Env_Temp", sensor="Humidity Sensor", channel_type="analog_in",
        sample_rate_hz=1, resolution_bits=12, voltage_range_v=(0, 3.3),
        notes="Ambient temperature"
    ))
    for axis in ['X', 'Y', 'Z']:
        channels.append(DAQChannel(
            name=f"Earth_B_{axis}", sensor="Earth Magnetometer", channel_type="analog_in",
            sample_rate_hz=100, resolution_bits=16, voltage_range_v=(-5, 5),
            notes=f"Earth field {axis}"
        ))
    
    # Digital I/O for control
    channels.append(DAQChannel(
        name="HV_Enable", sensor="Control", channel_type="digital_out",
        sample_rate_hz=1, resolution_bits=1, voltage_range_v=(0, 5),
        notes="HV supply enable"
    ))
    channels.append(DAQChannel(
        name="PWM_Drive", sensor="Control", channel_type="digital_out",
        sample_rate_hz=1000000, resolution_bits=1, voltage_range_v=(0, 5),
        notes="Gate drive PWM"
    ))
    channels.append(DAQChannel(
        name="E_Stop_Status", sensor="E-Stop", channel_type="digital_in",
        sample_rate_hz=1000, resolution_bits=1, voltage_range_v=(0, 5),
        notes="E-stop status monitor"
    ))
    channels.append(DAQChannel(
        name="Interlock_Status", sensor="Interlocks", channel_type="digital_in",
        sample_rate_hz=1000, resolution_bits=16, voltage_range_v=(0, 5),
        notes="Interlock status bitfield"
    ))
    
    return channels


def design_enclosure() -> Dict:
    """Design mechanical enclosure."""
    return {
        "type": "Ruggedized aluminum enclosure (IP65)",
        "dimensions_mm": {"width": 350, "depth": 280, "height": 200},
        "volume_l": 19.6,
        "mass_kg": 2.5,
        "material": "6061-T6 aluminum, anodized",
        "sealing": "Silicone gasket, IP65",
        "mounting": "4x M8 vibration-isolated mounts",
        "internal_layout": {
            "coil_chamber": "Central, 180mm dia x 60mm",
            "power_electronics": "Side compartment, forced air",
            "battery": "Bottom compartment, separate ventilation",
            "daq_compute": "Top compartment, conduction cooled"
        },
        "cable_glands": "8x M20 (power), 4x M12 (signal), 2x M16 (fiber)",
        "cooling_interface": "2x 80mm fan mounts + liquid cooling quick-connects",
        "cost_usd": 850
    }


def design_cooling() -> Dict:
    """Design cooling system."""
    return {
        "primary": "Forced air (2x 80mm PWM fans, 5 m/s)",
        "airflow_cfm": 45,
        "thermal_resistance_k_w": 0.20,
        "fan_power_w": 5,
        "fan_control": "PWM, temperature-regulated",
        "secondary": "Liquid cooling ready (quick-connects, 1/4\" tubing)",
        "liquid_flow_lpm": 2.0,
        "liquid_thermal_resistance_k_w": 0.02,
        "coolant": "50/50 ethylene glycol/water",
        "radiator": "External 120mm radiator (optional)",
        "temperature_monitoring": "12x RTD + 8x TC + IR camera",
        "overtemp_shutdown": "145 C winding, 55 C battery",
        "cost_usd": 180
    }


def design_test_procedures() -> List[Dict]:
    """Define test procedures."""
    return [
        {
            "id": "TP-001",
            "name": "Pre-Test Safety Verification",
            "description": "Verify all interlocks, E-stop, grounding, insulation resistance",
            "duration_min": 30,
            "prerequisites": ["Visual inspection", "Continuity test", "Hipot test (500V DC)"],
            "pass_criteria": "All interlocks functional, R_iso > 100 Mohm, ground < 0.1 ohm"
        },
        {
            "id": "TP-002",
            "name": "Sensor Calibration",
            "description": "Calibrate all sensors against traceable references",
            "duration_min": 120,
            "prerequisites": ["Calibration equipment", "Reference standards"],
            "pass_criteria": "All sensors within spec accuracy; calibration certs recorded"
        },
        {
            "id": "TP-003",
            "name": "Coil Characterization (Cold)",
            "description": "Measure R, L, Q at low current (1A) vs frequency",
            "duration_min": 45,
            "prerequisites": ["LCR meter", "Function generator"],
            "pass_criteria": "R_dc = 12.3 mOhm +/- 5%, L = 78.5 nH +/- 5%, Q > 50 at 1 MHz"
        },
        {
            "id": "TP-004",
            "name": "Power Supply Characterization",
            "description": "Verify HV supply output, PF correction, resonant operation",
            "duration_min": 60,
            "prerequisites": ["Power analyzer", "Oscilloscope", "Resistive load"],
            "pass_criteria": "V_out = 49.4V +/- 2%, PF > 0.95 with PFC, ZVS confirmed"
        },
        {
            "id": "TP-005",
            "name": "Low-Power Functional Test",
            "description": "Full system test at 10% current (10A), 1% duty",
            "duration_min": 30,
            "prerequisites": ["TP-001 to TP-004 passed"],
            "pass_criteria": "All sensors reading, interlocks active, waveforms clean"
        },
        {
            "id": "TP-006",
            "name": "Thermal Characterization",
            "description": "Steady-state thermal at 25%, 50%, 75% duty cycle",
            "duration_min": 240,
            "prerequisites": ["TP-005 passed", "Thermal equilibrium at each step"],
            "pass_criteria": "T_winding < 145C at 50% duty, thermal model validated"
        },
        {
            "id": "TP-007",
            "name": "Full-Power Pulse Test",
            "description": "Single 1ms pulse at 100A, capture all waveforms",
            "duration_min": 15,
            "prerequisites": ["TP-006 passed"],
            "pass_criteria": "I_peak = 100A +/- 2%, V_peak = 49.4V +/- 5%, no interlock trips"
        },
        {
            "id": "TP-008",
            "name": "Magnetic Field Mapping",
            "description": "3D B-field map at 100A, compare to simulation",
            "duration_min": 90,
            "prerequisites": ["TP-007 passed", "Hall array calibrated"],
            "pass_criteria": "B_center = 25 G +/- 5%, gradient = 250 G/m +/- 10%"
        },
        {
            "id": "TP-009",
            "name": "Extended Run Test",
            "description": "Continuous operation at 50% duty for 30 minutes",
            "duration_min": 45,
            "prerequisites": ["TP-008 passed"],
            "pass_criteria": "Stable temperatures, no drift, battery SoC tracked"
        },
        {
            "id": "TP-010",
            "name": "EMI/EMC Pre-Compliance",
            "description": "Radiated/conducted emissions scan",
            "duration_min": 120,
            "prerequisites": ["TP-009 passed", "EMI test facility"],
            "pass_criteria": "Below FCC Part 15 Class A limits with 6 dB margin"
        },
        {
            "id": "TP-011",
            "name": "Data Validation & Archival",
            "description": "Verify data integrity, checksum, archive raw + processed",
            "duration_min": 30,
            "prerequisites": ["All test data collected"],
            "pass_criteria": "Zero data loss, checksums match, metadata complete"
        },
    ]


def build_experimental_apparatus() -> ExperimentalApparatus:
    """Build complete experimental apparatus specification."""
    
    # Phase 10 coil spec (reference)
    coil_spec = {
        "wire_gauge": "AWG8",
        "major_radius_mm": 80,
        "minor_radius_mm": 10,
        "turns": 10,
        "current_a": 100,
        "frequency_hz": 1e6,
        "inductance_nh": 78.5,
        "resistance_mohm": 12.3,
        "mass_kg": 0.38,
        "volume_l": 0.16
    }
    
    hv_supply_spec = {
        "topology": "Class E resonant with PFC",
        "input_voltage_v": 48,
        "output_voltage_v": 49.4,
        "output_current_a": 100,
        "switching_freq_hz": 1e6,
        "efficiency": 0.92,
        "pfc_capacitor_nf": 323,
        "resonant_capacitor_nf": 323,
        "gaN_fets": "4x EPC2045 (100V, 35A)",
        "driver": "LMG1210 half-bridge driver",
        "control": "STM32H7 + FPGA (pulse shaping)",
        "mass_kg": 1.5,
        "volume_l": 1.0
    }
    
    battery_spec = {
        "chemistry": "Li-ion NMC 21700",
        "configuration": "13S4P (48V nominal, 20Ah)",
        "capacity_wh": 500,
        "max_continuous_a": 60,
        "peak_a": 120,
        "mass_kg": 3.5,
        "volume_l": 2.5,
        "bms": "Active balancing, CAN bus",
        "cycle_life": "1000 cycles to 80%",
        "cost_usd": 600
    }
    
    magnetic_sensors = design_magnetic_sensors()
    thermal_sensors = design_thermal_sensors()
    electrical_sensors = design_electrical_sensors()
    mechanical_sensors = design_mechanical_sensors()
    environmental_sensors = design_environmental_sensors()
    
    interlocks = design_safety_interlocks()
    emergency_stop = design_emergency_stop()
    
    daq_channels = design_daq_channels()
    
    enclosure = design_enclosure()
    cooling = design_cooling()
    
    test_procedures = design_test_procedures()
    
    # Calculate totals
    sensor_mass = sum(s.mass_g for s in magnetic_sensors + thermal_sensors + 
                      electrical_sensors + mechanical_sensors + environmental_sensors) / 1000
    sensor_power = sum(s.power_w for s in magnetic_sensors + thermal_sensors + 
                       electrical_sensors + mechanical_sensors + environmental_sensors)
    
    total_mass = (coil_spec["mass_kg"] + hv_supply_spec["mass_kg"] + 
                  battery_spec["mass_kg"] + enclosure["mass_kg"] + 
                  sensor_mass + 1.5)  # 1.5 kg structure/harness
    
    total_volume = (coil_spec["volume_l"] + hv_supply_spec["volume_l"] + 
                    battery_spec["volume_l"] + enclosure["volume_l"])
    
    # DAQ storage estimate
    total_sample_rate = sum(c.sample_rate_hz for c in daq_channels)
    # Assume 4 bytes per sample (32-bit float)
    storage_gb_per_hour = (total_sample_rate * 4 * 3600) / 1e9
    
    return ExperimentalApparatus(
        coil_spec=coil_spec,
        hv_supply_spec=hv_supply_spec,
        battery_spec=battery_spec,
        pfc_capacitor_f=323e-9,
        resonant_capacitor_f=323e-9,
        magnetic_sensors=magnetic_sensors,
        thermal_sensors=thermal_sensors,
        electrical_sensors=electrical_sensors,
        mechanical_sensors=mechanical_sensors,
        environmental_sensors=environmental_sensors,
        interlocks=interlocks,
        emergency_stop=emergency_stop,
        daq_channels=daq_channels,
        daq_sample_rate_hz=max(c.sample_rate_hz for c in daq_channels),
        daq_buffer_size=1000000,
        storage_gb_per_hour=storage_gb_per_hour,
        enclosure_spec=enclosure,
        cooling_spec=cooling,
        total_mass_kg=total_mass,
        total_volume_l=total_volume,
        test_procedures=test_procedures,
        safety_standards=[
            "IEC 61010-1 (Measurement equipment safety)",
            "IEC 60601-1 (Medical electrical equipment - reference for insulation)",
            "UL 61010-1",
            "NFPA 70E (Electrical safety in workplace)",
            "ISO 13849-1 (Safety-related control systems, PL d)"
        ],
        emc_standards=[
            "FCC Part 15 Class A (Radiated/conducted emissions)",
            "IEC 61000-4-2 (ESD immunity)",
            "IEC 61000-4-3 (Radiated immunity)",
            "IEC 61000-4-4 (EFT/Burst immunity)",
            "IEC 61000-4-5 (Surge immunity)",
            "IEC 61000-4-6 (Conducted immunity)",
            "MIL-STD-461G (Military EMC - reference)"
        ]
    )


def save_apparatus_spec(apparatus: ExperimentalApparatus, filepath: str):
    """Save apparatus specification to JSON."""
    # Convert to serializable dict
    def dc_to_dict(obj):
        if hasattr(obj, '__dataclass_fields__'):
            return {k: dc_to_dict(v) for k, v in obj.__dict__.items()}
        elif isinstance(obj, list):
            return [dc_to_dict(v) for v in obj]
        elif isinstance(obj, dict):
            return {k: dc_to_dict(v) for k, v in obj.items()}
        elif isinstance(obj, tuple):
            return list(obj)
        else:
            return obj
    
    data = dc_to_dict(apparatus)
    
    with open(filepath, 'w') as f:
        json.dump(data, f, indent=2)
    
    print(f"Saved apparatus spec to {filepath}")


def print_summary(apparatus: ExperimentalApparatus):
    """Print apparatus summary."""
    print("=" * 70)
    print("PHASE 11: EXPERIMENTAL APPARATUS DESIGN SUMMARY")
    print("=" * 70)
    
    print(f"\n--- SYSTEM TOTALS ---")
    print(f"Total Mass: {apparatus.total_mass_kg:.1f} kg (budget: 10 kg)")
    print(f"Total Volume: {apparatus.total_volume_l:.1f} L (budget: 20 L)")
    print(f"Mass Margin: {10 - apparatus.total_mass_kg:.1f} kg")
    print(f"Volume Margin: {20 - apparatus.total_volume_l:.1f} L")
    
    print(f"\n--- COIL SYSTEM (Phase 10) ---")
    cs = apparatus.coil_spec
    print(f"  AWG{cs['wire_gauge']}, R={cs['major_radius_mm']}mm, a={cs['minor_radius_mm']}mm, N={cs['turns']}")
    print(f"  {cs['current_a']}A @ {cs['frequency_hz']/1e6:.0f}MHz, L={cs['inductance_nh']:.1f}nH")
    
    print(f"\n--- POWER SYSTEM ---")
    hs = apparatus.hv_supply_spec
    print(f"  Topology: {hs['topology']}")
    print(f"  {hs['input_voltage_v']}V -> {hs['output_voltage_v']}V @ {hs['output_current_a']}A")
    print(f"  PFC Cap: {apparatus.pfc_capacitor_f*1e9:.0f} nF, Resonant Cap: {apparatus.resonant_capacitor_f*1e9:.0f} nF")
    
    bs = apparatus.battery_spec
    print(f"  Battery: {bs['chemistry']} {bs['configuration']}, {bs['capacity_wh']} Wh")
    print(f"  Mass: {bs['mass_kg']} kg, Volume: {bs['volume_l']} L")
    
    print(f"\n--- SENSORS ---")
    all_sensors = (apparatus.magnetic_sensors + apparatus.thermal_sensors + 
                   apparatus.electrical_sensors + apparatus.mechanical_sensors + 
                   apparatus.environmental_sensors)
    print(f"  Total sensors: {len(all_sensors)}")
    print(f"  Magnetic: {len(apparatus.magnetic_sensors)}")
    print(f"  Thermal: {len(apparatus.thermal_sensors)}")
    print(f"  Electrical: {len(apparatus.electrical_sensors)}")
    print(f"  Mechanical: {len(apparatus.mechanical_sensors)}")
    print(f"  Environmental: {len(apparatus.environmental_sensors)}")
    print(f"  Sensor mass: {sum(s.mass_g for s in all_sensors)/1000:.2f} kg")
    print(f"  Sensor power: {sum(s.power_w for s in all_sensors):.2f} W")
    
    print(f"\n--- SAFETY INTERLOCKS ---")
    print(f"  Count: {len(apparatus.interlocks)} + 1 E-Stop")
    for il in apparatus.interlocks:
        print(f"  - {il.name}: {il.trigger_condition} -> {il.action} (RL{il.redundancy_level})")
    print(f"  - E-Stop: {apparatus.emergency_stop.response_time_ms}ms hard cutoff")
    
    print(f"\n--- DATA ACQUISITION ---")
    print(f"  Channels: {len(apparatus.daq_channels)}")
    print(f"  Max sample rate: {apparatus.daq_sample_rate_hz/1e6:.1f} MHz")
    print(f"  Storage: {apparatus.storage_gb_per_hour:.2f} GB/hour")
    print(f"  Buffer: {apparatus.daq_buffer_size} samples")
    
    print(f"\n--- MECHANICAL ---")
    enc = apparatus.enclosure_spec
    print(f"  Enclosure: {enc['type']}")
    print(f"  Dimensions: {enc['dimensions_mm']['width']}x{enc['dimensions_mm']['depth']}x{enc['dimensions_mm']['height']} mm")
    print(f"  Volume: {enc['volume_l']:.1f} L, Mass: {enc['mass_kg']} kg")
    
    cool = apparatus.cooling_spec
    print(f"  Cooling: {cool['primary']}")
    print(f"  Thermal R: {cool['thermal_resistance_k_w']:.2f} K/W")
    print(f"  Liquid ready: {cool['secondary']}")
    
    print(f"\n--- TEST PROCEDURES ---")
    for tp in apparatus.test_procedures:
        print(f"  {tp['id']}: {tp['name']} ({tp['duration_min']} min)")
    
    print(f"\n--- COMPLIANCE ---")
    print(f"  Safety: {', '.join(apparatus.safety_standards[:3])}...")
    print(f"  EMC: {', '.join(apparatus.emc_standards[:3])}...")


def main():
    """Main entry point."""
    print("Building Phase 11 Experimental Apparatus Specification...")
    
    apparatus = build_experimental_apparatus()
    
    print_summary(apparatus)
    
    # Save to results
    import os
    os.makedirs("results/validated", exist_ok=True)
    save_apparatus_spec(apparatus, "results/validated/phase11_experimental_apparatus.json")
    
    print("\n" + "=" * 70)
    print("Phase 11 Experimental Apparatus Design COMPLETE")
    print("=" * 70)
    
    return apparatus


if __name__ == "__main__":
    main()