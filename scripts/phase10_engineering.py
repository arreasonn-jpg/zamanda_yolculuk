#!/usr/bin/env python3
"""
Phase 10: Engineering Feasibility Analysis

Coil design optimization, thermal analysis, HV/power supply,
and backpack constraint optimizer for the time-travel device.
"""

import numpy as np
import json
from pathlib import Path
from dataclasses import dataclass, asdict
from typing import Dict, List, Tuple
import sys

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from zt005.model import Backpack, Drive
from zt005.geometry import build_sources


# Physical constants
MU_0 = 4 * np.pi * 1e-7  # H/m
COPPER_RESISTIVITY = 1.68e-8  # ohm*m at 20C
COPPER_DENSITY = 8960  # kg/m^3
COPPER_SPECIFIC_HEAT = 385  # J/(kg*K)
COPPER_TEMP_COEFF = 0.00393  # per K
INSULATION_THICKNESS = 0.1e-3  # 0.1 mm enamel
MAX_WIRE_TEMP = 150  # C (enamel limit)
AMBIENT_TEMP = 25  # C
MAX_SURFACE_TEMP = 60  # C (touch safety)


@dataclass
class CoilSpec:
    """Coil specification."""
    name: str
    # Geometry
    major_radius_m: float  # R (distance from center to wire center)
    minor_radius_m: float  # a (wire bundle radius)
    turns: int
    # Wire
    wire_diameter_mm: float  # bare copper diameter
    insulation: str = "enamel"
    # Electrical
    current_a: float = 100.0
    frequency_hz: float = 1e6
    # Derived
    length_m: float = 0.0
    resistance_ohm: float = 0.0
    inductance_h: float = 0.0
    mass_kg: float = 0.0
    power_w: float = 0.0
    voltage_v: float = 0.0
    # Thermal
    temp_rise_c: float = 0.0
    duty_cycle_max: float = 1.0
    # Field
    b_field_t: float = 0.0
    b_field_gradient_t_m: float = 0.0


def compute_coil_spec(spec: CoilSpec) -> CoilSpec:
    """Compute derived coil parameters."""
    # Wire cross-section (bare copper)
    d = spec.wire_diameter_mm * 1e-3
    a_wire = np.pi * (d/2)**2
    
    # Total wire length (approximate as circular loops)
    spec.length_m = spec.turns * 2 * np.pi * spec.major_radius_m
    
    # Resistance at operating temperature (estimate 75C average)
    rho_75 = COPPER_RESISTIVITY * (1 + COPPER_TEMP_COEFF * (75 - 20))
    spec.resistance_ohm = rho_75 * spec.length_m / a_wire
    
    # Inductance (single layer solenoid approximation)
    # L = mu0 * N^2 * A / l_eff where l_eff ~ 2*pi*R for toroid
    # For a toroidal coil: L = mu0 * N^2 * (minor_radius^2) / (2 * major_radius)
    spec.inductance_h = MU_0 * spec.turns**2 * spec.minor_radius_m**2 / (2 * spec.major_radius_m)
    
    # Copper mass
    spec.mass_kg = COPPER_DENSITY * a_wire * spec.length_m
    
    # Power (I^2 * R)
    spec.power_w = spec.current_a**2 * spec.resistance_ohm
    
    # Voltage (V = I * R + I * omega * L for AC)
    omega = 2 * np.pi * spec.frequency_hz
    spec.voltage_v = spec.current_a * np.sqrt(spec.resistance_ohm**2 + (omega * spec.inductance_h)**2)
    
    # Magnetic field at center of toroid
    # B = mu0 * N * I / (2 * pi * R)
    spec.b_field_t = MU_0 * spec.turns * spec.current_a / (2 * np.pi * spec.major_radius_m)
    
    # Field gradient (approximate)
    spec.b_field_gradient_t_m = spec.b_field_t / spec.minor_radius_m
    
    # Thermal: steady-state temperature rise
    # Surface area for cooling (cylinder approximation)
    surface_area = 2 * np.pi * spec.major_radius_m * (2 * spec.minor_radius_m)  # rough
    # Convection coefficient (natural convection in air) ~ 10 W/(m^2*K)
    h_conv = 10.0
    spec.temp_rise_c = spec.power_w / (h_conv * surface_area)
    
    # Max duty cycle for temperature limit
    max_temp_rise = MAX_WIRE_TEMP - AMBIENT_TEMP
    if spec.temp_rise_c > 0:
        spec.duty_cycle_max = min(1.0, max_temp_rise / spec.temp_rise_c)
    else:
        spec.duty_cycle_max = 1.0
    
    return spec


def optimize_coil_for_backpack(
    max_weight_kg: float = 10.0,
    max_volume_l: float = 20.0,
    max_power_w: float = 200.0,
    target_b_field_t: float = 1.0,
    current_a: float = 100.0,
    frequency_hz: float = 1e6
) -> List[CoilSpec]:
    """
    Optimize coil design for backpack constraints.
    
    Returns list of viable designs sorted by figure of merit.
    """
    results = []
    
    # Wire gauge options (AWG) - use thicker wire for lower resistance
    awg_options = [
        (8, 3.264), (10, 2.588), (12, 2.053), (14, 1.628), (16, 1.291),
        (18, 1.024), (20, 0.812), (22, 0.644), (24, 0.511)
    ]
    
    # Major radius range (coil size) - larger radius reduces inductance
    major_radii = np.linspace(0.08, 0.25, 18)  # 8-25 cm
    
    # Minor radius (coil thickness)
    minor_radii = np.linspace(0.01, 0.05, 9)  # 1-5 cm
    
    # Turns range - fewer turns for lower inductance at 1MHz
    turns_range = range(10, 201, 10)
    
    count = 0
    for awg, dia_mm in awg_options:
        for R in major_radii:
            for a in minor_radii:
                for N in turns_range:
                    spec = CoilSpec(
                        name=f"AWG{awg}_R{R*100:.0f}mm_a{a*100:.0f}mm_N{N}",
                        major_radius_m=R,
                        minor_radius_m=a,
                        turns=N,
                        wire_diameter_mm=dia_mm,
                        current_a=current_a,
                        frequency_hz=frequency_hz
                    )
                    spec = compute_coil_spec(spec)
                    
                    # Check constraints
                    # Volume (approximate as torus)
                    volume_m3 = (2 * np.pi * R) * (np.pi * a**2)
                    volume_l = volume_m3 * 1000
                    
                    # Debug first few
                    if count < 5:
                        print(f"  Debug: {spec.name}: mass={spec.mass_kg:.3f}kg, vol={volume_l:.3f}L, power={spec.power_w:.1f}W, duty={spec.duty_cycle_max:.3f}, V={spec.voltage_v:.1f}V")
                    
                    if spec.mass_kg > max_weight_kg:
                        continue
                    if volume_l > max_volume_l:
                        continue
                    if spec.power_w > max_power_w:
                        continue
                    if spec.duty_cycle_max < 0.1:  # At least 10% duty cycle
                        continue
                    if spec.voltage_v > 1000:  # Safety limit
                        continue
                    
                    # Figure of merit: B-field per sqrt(power) per weight
                    fom = spec.b_field_t / np.sqrt(spec.power_w * spec.mass_kg)
                    
                    results.append((fom, spec))
                    count += 1
    
    print(f"  Total combinations checked: {count}, viable designs: {len(results)}")
    
    # Sort by figure of merit (descending)
    results.sort(key=lambda x: x[0], reverse=True)
    
    return [spec for fom, spec in results]


def thermal_analysis(spec: CoilSpec, duty_cycle: float = 1.0, 
                     cooling: str = "natural") -> Dict:
    """
    Detailed thermal analysis.
    
    cooling: "natural", "forced_air", "liquid"
    """
    # Heat generation
    p_avg = spec.power_w * duty_cycle
    
    # Thermal resistance
    surface_area = 2 * np.pi * spec.major_radius_m * (2 * spec.minor_radius_m)
    
    if cooling == "natural":
        h = 10.0  # W/(m^2*K)
    elif cooling == "forced_air":
        h = 50.0  # W/(m^2*K) with fan
    elif cooling == "liquid":
        h = 500.0  # W/(m^2*K) with liquid cooling
    else:
        h = 10.0
    
    r_th = 1.0 / (h * surface_area)  # K/W
    
    # Steady state
    temp_rise_steady = p_avg * r_th
    temp_steady = AMBIENT_TEMP + temp_rise_steady
    
    # Thermal time constant
    thermal_mass = spec.mass_kg * COPPER_SPECIFIC_HEAT  # J/K
    tau = thermal_mass * r_th  # seconds
    
    # Transient (step response)
    def temp_at_time(t):
        return AMBIENT_TEMP + temp_rise_steady * (1 - np.exp(-t / tau))
    
    return {
        "power_avg_w": p_avg,
        "thermal_resistance_k_w": r_th,
        "thermal_time_constant_s": tau,
        "steady_state_temp_c": temp_steady,
        "temp_rise_steady_c": temp_rise_steady,
        "temp_at_10s_c": temp_at_time(10),
        "temp_at_60s_c": temp_at_time(60),
        "temp_at_300s_c": temp_at_time(300),
        "max_safe_duty_cycle": min(1.0, (MAX_WIRE_TEMP - AMBIENT_TEMP) / temp_rise_steady) if temp_rise_steady > 0 else 1.0,
        "surface_temp_c": AMBIENT_TEMP + p_avg / (h * surface_area) * 0.5  # approximate
    }


def hv_power_supply_design(spec: CoilSpec) -> Dict:
    """Design HV power supply specifications."""
    omega = 2 * np.pi * spec.frequency_hz
    
    # Impedance
    z_mag = np.sqrt(spec.resistance_ohm**2 + (omega * spec.inductance_h)**2)
    phase_deg = np.degrees(np.arctan2(omega * spec.inductance_h, spec.resistance_ohm))
    
    # Reactive power
    q_var = spec.current_a**2 * omega * spec.inductance_h
    
    # Apparent power
    s_va = spec.current_a * spec.voltage_v
    
    # Power factor
    pf = spec.power_w / s_va if s_va > 0 else 0
    
    # Capacitor for power factor correction (parallel)
    c_pf = 1 / (omega * z_mag * np.sin(np.radians(phase_deg))) if phase_deg > 0 else 0
    
    # Resonant capacitor (series resonance)
    c_res = 1 / (omega**2 * spec.inductance_h)
    
    return {
        "impedance_ohm": z_mag,
        "phase_deg": phase_deg,
        "real_power_w": spec.power_w,
        "reactive_power_var": q_var,
        "apparent_power_va": s_va,
        "power_factor": pf,
        "voltage_v": spec.voltage_v,
        "current_a": spec.current_a,
        "pf_correction_capacitor_f": c_pf,
        "series_resonant_capacitor_f": c_res,
        "insulation_class": "Class F (155C)" if spec.voltage_v < 1000 else "Class H (180C)",
        "creepage_distance_mm": max(2.5, spec.voltage_v / 400),  # IEC 60664 approx
        "clearance_mm": max(1.5, spec.voltage_v / 600)
    }


def backpack_optimizer(
    configs: List[str] = ["G1", "G2", "G3", "G4"],
    max_weight_kg: float = 10.0,
    max_volume_l: float = 20.0,
    max_power_w: float = 200.0,
    battery_wh: float = 500.0,
    current_a: float = 100.0,
    frequency_hz: float = 1e6
) -> Dict:
    """
    Full backpack system optimization.
    
    Returns best design for each config.
    """
    backpack = Backpack()
    results = {}
    
    for config_name in configs:
        sources = build_sources(config_name, backpack)
        
        # Get coil parameters from source geometry
        # For simplicity, use the first loop as reference
        import numpy as np
        if sources:
            pts, seg = sources[0]  # sources is a list of (pts, seg) tuples
            # Compute radius from points (distance from center)
            center = np.mean(pts, axis=0)
            radii = np.linalg.norm(pts - center, axis=1)
            R = np.mean(radii)  # major radius
            # Estimate minor radius from wire bundle
            a = 0.02  # 2 cm default
        else:
            R = 0.1
            a = 0.02
        
        # Optimize
        designs = optimize_coil_for_backpack(
            max_weight_kg=max_weight_kg,
            max_volume_l=max_volume_l,
            max_power_w=max_power_w,
            current_a=current_a,
            frequency_hz=frequency_hz
        )
        
        if designs:
            best = designs[0]
            thermal = thermal_analysis(best, duty_cycle=0.5, cooling="forced_air")
            hv = hv_power_supply_design(best)
            
            # Battery life
            energy_per_pulse_j = best.power_w * 1.0  # 1 second pulse
            battery_j = battery_wh * 3600
            max_pulses = battery_j / energy_per_pulse_j if energy_per_pulse_j > 0 else 0
            
            results[config_name] = {
                "coil": asdict(best),
                "thermal": thermal,
                "hv_supply": hv,
                "battery": {
                    "capacity_wh": battery_wh,
                    "energy_per_pulse_j": energy_per_pulse_j,
                    "max_pulses": max_pulses,
                    "runtime_at_50_duty_s": battery_j / (best.power_w * 0.5) if best.power_w > 0 else 0
                },
                "alternatives": [asdict(d) for d in designs[1:4]]  # Top 3 alternatives
            }
        else:
            results[config_name] = {"error": "No viable design found"}
    
    return results


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Phase 10: Engineering Feasibility")
    parser.add_argument("--configs", nargs="+", default=["G1", "G2", "G3", "G4"])
    parser.add_argument("--max-weight", type=float, default=10.0)
    parser.add_argument("--max-volume", type=float, default=20.0)
    parser.add_argument("--max-power", type=float, default=200.0)
    parser.add_argument("--battery-wh", type=float, default=500.0)
    parser.add_argument("--current", type=float, default=100.0)
    parser.add_argument("--freq", type=float, default=1e6)
    parser.add_argument("--output-dir", default="results/validated")
    args = parser.parse_args()
    
    print("=" * 70)
    print("PHASE 10: ENGINEERING FEASIBILITY ANALYSIS")
    print("=" * 70)
    print(f"Configs: {args.configs}")
    print(f"Constraints: {args.max_weight} kg, {args.max_volume} L, {args.max_power} W")
    print(f"Battery: {args.battery_wh} Wh")
    print(f"Drive: {args.current} A @ {args.freq/1e6:.1f} MHz")
    print()
    
    results = backpack_optimizer(
        configs=args.configs,
        max_weight_kg=args.max_weight,
        max_volume_l=args.max_volume,
        max_power_w=args.max_power,
        battery_wh=args.battery_wh,
        current_a=args.current,
        frequency_hz=args.freq
    )
    
    # Print summary
    for config_name, data in results.items():
        if "error" in data:
            print(f"\n{config_name}: {data['error']}")
            continue
        
        coil = data["coil"]
        thermal = data["thermal"]
        hv = data["hv_supply"]
        battery = data["battery"]
        
        print(f"\n{'='*70}")
        print(f"CONFIG: {config_name}")
        print(f"{'='*70}")
        print(f"Coil: {coil['name']}")
        print(f"  Geometry: R={coil['major_radius_m']*100:.1f}mm, a={coil['minor_radius_m']*100:.1f}mm, N={coil['turns']}")
        print(f"  Wire: AWG ({coil['wire_diameter_mm']:.2f}mm), L={coil['length_m']:.1f}m, m={coil['mass_kg']:.2f}kg")
        print(f"  Electrical: R={coil['resistance_ohm']:.3f} ohm, L={coil['inductance_h']*1e6:.1f} uH")
        print(f"  Power: {coil['power_w']:.1f}W, V={coil['voltage_v']:.1f}V (RMS)")
        print(f"  Field: B={coil['b_field_t']*1e4:.2f}G, gradB={coil['b_field_gradient_t_m']:.1f}T/m")
        print(f"  Thermal: dT={coil['temp_rise_c']:.1f} C, duty_max={coil['duty_cycle_max']*100:.0f}%")
        print(f"  Thermal (forced air, 50% duty): steady={thermal['steady_state_temp_c']:.1f} C, tau={thermal['thermal_time_constant_s']:.1f}s")
        print(f"  HV Supply: Z={hv['impedance_ohm']:.2f} ohm, PF={hv['power_factor']:.3f}, S={hv['apparent_power_va']:.0f}VA")
        print(f"  Battery: {battery['max_pulses']:.0f} pulses, {battery['runtime_at_50_duty_s']/60:.1f}min @ 50% duty")
    
    # Save results
    output_path = Path(args.output_dir) / "phase10_engineering_feasibility.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults written to {output_path}")
    
    return results


if __name__ == "__main__":
    main()