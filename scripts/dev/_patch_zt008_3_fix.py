"""Fix ZT-008.3: add hardware/mechanical layers to F05, F06."""
from pathlib import Path

SRC = Path("src/zt005/run_zt008_3.py")
c = SRC.read_text(encoding="utf-8")

# F05: Communication loss -> add hardware dead-man relay
old_f05 = '''        {
            "id": "F05",
            "component": "Communication loss (BLE/WiFi)",
            "cause": "Radio interference / app crash",
            "severity": "low",
            "detection": "Timeout > 3 s",
            "mitigation": "Auto-shutdown device",
            "layer": "firmware",
        },'''

new_f05 = '''        {
            "id": "F05",
            "component": "Communication loss (BLE/WiFi)",
            "cause": "Radio interference / app crash",
            "severity": "low",
            "detection": "Timeout > 3 s (firmware) OR dead-man relay open",
            "mitigation": "Firmware auto-shutdown + hardware dead-man relay "
                          "(requires periodic heartbeat pulse)",
            "layer": "firmware + hardware",
        },'''

# F06: Mechanical breach -> add mechanical shock switch
old_f06 = '''        {
            "id": "F06",
            "component": "Mechanical enclosure breach",
            "cause": "Drop / impact",
            "severity": "medium",
            "detection": "IMU spike > 5g",
            "mitigation": "Auto-shutdown + safety review",
            "layer": "firmware",
        },'''

new_f06 = '''        {
            "id": "F06",
            "component": "Mechanical enclosure breach",
            "cause": "Drop / impact",
            "severity": "medium",
            "detection": "IMU spike > 5g (firmware) OR mechanical shock "
                         "switch trips at > 10g",
            "mitigation": "Firmware auto-shutdown + mechanical shock switch "
                          "hardwired to coil relay",
            "layer": "firmware + mechanical",
        },'''

assert old_f05 in c, "F05 block not found"
assert old_f06 in c, "F06 block not found"
c = c.replace(old_f05, new_f05).replace(old_f06, new_f06)
SRC.write_text(c, encoding="utf-8")
print("[OK] Fixed F05 (hardware) and F06 (mechanical)")
print()
print("Now run:")
print("  python -m pytest tests\\test_zt008_3.py -v")
print("  python -m zt005.run_zt008_3 --mode quick")