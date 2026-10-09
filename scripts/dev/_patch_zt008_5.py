"""Auto-patch for ZT-008.5 — Mobile app specification."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src" / "zt005" / "run_zt008_5.py"
TST = ROOT / "tests" / "test_zt008_5.py"

SRC_CODE = r'''
"""ZT-008.5 - Mobile application specification."""
import argparse, json
from pathlib import Path
import numpy as np


def tech_stack_options():
    """Candidate mobile frameworks with pros/cons."""
    return [
        {
            "name": "Flutter (Dart)",
            "platforms": ["Android", "iOS"],
            "ble": "flutter_blue_plus",
            "wifi": "http/websocket",
            "usb": "usb_serial",
            "score": 9,
            "pros": ["Single codebase", "Mature BLE", "Good UI perf"],
            "cons": ["Dart learning curve"],
            "recommendation": "PRIMARY",
        },
        {
            "name": "React Native (JS)",
            "platforms": ["Android", "iOS"],
            "ble": "react-native-ble-plx",
            "wifi": "fetch/websocket",
            "usb": "react-native-usb-serial",
            "score": 7,
            "pros": ["JS ecosystem", "Fast dev"],
            "cons": ["Bridge performance issues", "USB immaturity"],
            "recommendation": "FALLBACK",
        },
        {
            "name": "Native (Kotlin + Swift)",
            "platforms": ["Android", "iOS"],
            "ble": "native APIs",
            "wifi": "native",
            "usb": "native",
            "score": 10,
            "pros": ["Best perf", "Full API access"],
            "cons": ["2x codebase", "Slow dev"],
            "recommendation": "NO (too slow)",
        },
        {
            "name": "Progressive Web App (PWA)",
            "platforms": ["Web"],
            "ble": "Web Bluetooth (limited)",
            "wifi": "native",
            "usb": "WebUSB (limited)",
            "score": 5,
            "pros": ["No install"],
            "cons": ["Limited BLE", "No USB on iOS"],
            "recommendation": "NO (insufficient)",
        },
    ]


def screen_list():
    """Screens and their purpose."""
    return [
        {
            "name": "Splash",
            "purpose": "Init + device discovery",
            "state_data": [],
            "actions": ["scan", "connect"],
        },
        {
            "name": "Home / Dashboard",
            "purpose": "Live telemetry view",
            "state_data": ["link_status", "battery", "power",
                           "temperature", "coil_current"],
            "actions": ["arm", "view_logs", "settings"],
        },
        {
            "name": "Target Selector",
            "purpose": "Date/time/location input for temporal target",
            "state_data": ["target_date", "target_time",
                           "target_lat", "target_lon", "delta_t"],
            "actions": ["validate", "save_target"],
        },
        {
            "name": "Simulation Preview",
            "purpose": "Show predicted metric effect + safety state",
            "state_data": ["h_00_predicted", "delta_tau_per_year",
                           "energy_required", "safety_flags"],
            "actions": ["proceed_to_arm"],
        },
        {
            "name": "Arm / Ready",
            "purpose": "Two-factor authorization",
            "state_data": ["physical_key_present", "app_auth_ok",
                           "interlock_status"],
            "actions": ["confirm_arm", "abort"],
        },
        {
            "name": "Run Monitor",
            "purpose": "Live monitoring during operation",
            "state_data": ["elapsed_time", "all_sensors",
                           "fault_flags", "log_status"],
            "actions": ["emergency_stop", "stop_normal"],
        },
        {
            "name": "Post-Run Report",
            "purpose": "Summary + data export",
            "state_data": ["run_id", "duration", "log_file",
                           "anomalies"],
            "actions": ["export_csv", "share", "upload_to_lab"],
        },
        {
            "name": "Logs / History",
            "purpose": "Browse past runs",
            "state_data": ["run_list", "run_metadata"],
            "actions": ["view", "delete", "export"],
        },
        {
            "name": "Settings",
            "purpose": "Device + app configuration",
            "state_data": ["link_pref", "sample_rate", "clock_source",
                           "safety_limits"],
            "actions": ["save", "reset"],
        },
        {
            "name": "Help / Safety",
            "purpose": "Instructions + emergency info",
            "state_data": [],
            "actions": ["view_docs"],
        },
    ]


def ble_protocol():
    """BLE GATT service/characteristic spec."""
    return {
        "service_uuid": "b7e0_0001_zt_svc",  # placeholder
        "characteristics": [
            {
                "name": "Status",
                "uuid": "b7e0_0002",
                "access": "read + notify",
                "size_bytes": 16,
                "fields": ["state", "battery_pct", "temperature_C",
                           "flags"],
            },
            {
                "name": "Command",
                "uuid": "b7e0_0003",
                "access": "write",
                "size_bytes": 8,
                "commands": ["ARM", "START", "STOP", "ABORT", "RESET",
                             "SET_TARGET", "GET_STATUS"],
            },
            {
                "name": "Telemetry_Stream",
                "uuid": "b7e0_0004",
                "access": "notify",
                "size_bytes": 244,
                "note": "BLE 5.0 max payload; sent at 1 Hz for summary",
            },
            {
                "name": "Safety_Interlock",
                "uuid": "b7e0_0005",
                "access": "read",
                "size_bytes": 8,
                "fields": ["hw_estop", "watchdog_ok", "relay_closed",
                           "thermal_ok"],
            },
        ],
    }


def wifi_protocol():
    """WiFi REST + WebSocket spec."""
    return {
        "rest_endpoints": [
            {"method": "GET",  "path": "/api/status",
             "returns": "device state JSON"},
            {"method": "GET",  "path": "/api/sensors",
             "returns": "latest sensor values"},
            {"method": "POST", "path": "/api/arm",
             "returns": "arm ack"},
            {"method": "POST", "path": "/api/start",
             "returns": "run_id"},
            {"method": "POST", "path": "/api/stop",
             "returns": "final summary"},
            {"method": "POST", "path": "/api/abort",
             "returns": "immediate ack"},
            {"method": "GET",  "path": "/api/logs",
             "returns": "list of past runs"},
            {"method": "GET",  "path": "/api/logs/{run_id}",
             "returns": "log metadata"},
            {"method": "GET",  "path": "/api/logs/{run_id}/csv",
             "returns": "full CSV stream"},
        ],
        "websocket": {
            "endpoint": "ws://<device-ip>:8080/telemetry",
            "rate_hz": 100,
            "format": "protobuf or msgpack",
            "note": "Raw sensor stream during RUN",
        },
    }


def app_state_machine():
    """App-side state machine."""
    return {
        "DISCONNECTED": {"CONNECT_BLE": "DISCOVERING",
                         "CONNECT_WIFI": "DISCOVERING"},
        "DISCOVERING":  {"DEVICE_FOUND": "CONNECTED",
                         "TIMEOUT": "DISCONNECTED"},
        "CONNECTED":    {"ARM_REQ": "ARMING",
                         "DISCONNECT": "DISCONNECTED"},
        "ARMING":       {"AUTH_OK": "ARMED",
                         "AUTH_FAIL": "CONNECTED"},
        "ARMED":        {"START": "RUNNING",
                         "ABORT": "CONNECTED"},
        "RUNNING":      {"STOP": "POST_RUN",
                         "EMERGENCY_STOP": "POST_RUN",
                         "FAULT": "POST_RUN",
                         "LINK_LOSS": "FAULT_HANDLING"},
        "FAULT_HANDLING": {"RESYNC": "CONNECTED",
                            "ABORT": "POST_RUN"},
        "POST_RUN":     {"ACK": "CONNECTED"},
    }


def compute_link_requirements():
    """What does each screen need from the link?"""
    return {
        "Splash":         {"ble": True, "wifi": False, "usb": False},
        "Dashboard":      {"ble": True, "wifi": True,  "usb": False},
        "Run_Monitor":    {"ble": True, "wifi": True,  "usb": True},
        "Post_Run":       {"ble": False, "wifi": True, "usb": True},
        "Logs":           {"ble": False, "wifi": True, "usb": True},
    }


def save_screens_figure(out_path):
    """Draw screen flow diagram."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        screens = [s["name"] for s in screen_list()]
        n = len(screens)
        # Simple vertical flow
        fig, ax = plt.subplots(figsize=(8, 12))
        for i, name in enumerate(screens):
            y = n - i
            ax.add_patch(plt.Rectangle((0.2, y - 0.4), 0.6, 0.8,
                                        facecolor="#e6f3ff",
                                        edgecolor="#1f77b4"))
            ax.text(0.5, y, name, ha="center", va="center", fontsize=11)
            if i < n - 1:
                ax.annotate("", xy=(0.5, y - 0.9), xytext=(0.5, y - 0.5),
                            arrowprops=dict(arrowstyle="->",
                                            color="gray"))
        ax.set_xlim(0, 1)
        ax.set_ylim(0, n + 1)
        ax.axis("off")
        ax.set_title("Mobile app screen flow")
        plt.tight_layout()
        fig.savefig(out_path, dpi=120)
        plt.close(fig)
        return True
    except Exception as e:
        print(f"Figure save failed: {e}")
        return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="quick", choices=["quick", "full"])
    a = ap.parse_args()

    print("=" * 96)
    print("ZAMANDA YOLCULUK - ZT-008.5")
    print("MOBILE APP SPECIFICATION")
    print("=" * 96)
    print()

    # --- Tech stack ---
    print("=" * 96)
    print("TECH STACK CANDIDATES")
    print("=" * 96)
    print(f"{'Framework':<28s} {'Platforms':<20s} {'Score':>6s} "
          f"{'Recommendation':<16s}")
    print("-" * 96)
    for fw in tech_stack_options():
        print(f"{fw['name']:<28s} {','.join(fw['platforms']):<20s} "
              f"{fw['score']:>6d} {fw['recommendation']:<16s}")
    print()
    print("SECIM: Flutter (single codebase + mature BLE)")
    print()

    # --- Screens ---
    print("=" * 96)
    print("SCREEN LIST")
    print("=" * 96)
    for i, s in enumerate(screen_list(), 1):
        print(f"  {i:>2d}. {s['name']:<22s}  {s['purpose']}")
    print()

    # --- BLE protocol ---
    print("=" * 96)
    print("BLE GATT PROTOCOL")
    print("=" * 96)
    ble = ble_protocol()
    print(f"Service UUID: {ble['service_uuid']}")
    print()
    for ch in ble["characteristics"]:
        print(f"  {ch['name']:<20s} UUID: {ch['uuid']:<20s} "
              f"access={ch['access']:<15s} size={ch['size_bytes']}B")
        if "commands" in ch:
            print(f"    commands: {ch['commands']}")
        if "fields" in ch:
            print(f"    fields: {ch['fields']}")
    print()

    # --- WiFi protocol ---
    print("=" * 96)
    print("WiFi REST + WebSocket PROTOCOL")
    print("=" * 96)
    wifi = wifi_protocol()
    print("REST endpoints:")
    for ep in wifi["rest_endpoints"]:
        print(f"  {ep['method']:<6s} {ep['path']:<30s} -> {ep['returns']}")
    print()
    ws = wifi["websocket"]
    print(f"WebSocket: {ws['endpoint']}  ({ws['rate_hz']} Hz, "
          f"{ws['format']})")
    print()

    # --- App state machine ---
    print("=" * 96)
    print("APP STATE MACHINE")
    print("=" * 96)
    for state, trans in app_state_machine().items():
        print(f"  {state:<18s} -> {dict(trans)}")
    print()

    # --- Link requirements ---
    print("=" * 96)
    print("LINK REQUIREMENTS PER SCREEN")
    print("=" * 96)
    print(f"{'Screen':<20s} {'BLE':>6s} {'WiFi':>6s} {'USB':>6s}")
    print("-" * 40)
    for scr, req in compute_link_requirements().items():
        print(f"{scr:<20s} {str(req['ble']):>6s} "
              f"{str(req['wifi']):>6s} {str(req['usb']):>6s}")
    print()

    # --- Summary ---
    print("=" * 96)
    print("SONUC")
    print("=" * 96)
    print(f"  Recommended stack:   Flutter (Dart)")
    print(f"  Screens:             {len(screen_list())}")
    print(f"  BLE characteristics: {len(ble['characteristics'])}")
    print(f"  REST endpoints:      {len(wifi['rest_endpoints'])}")
    print(f"  App states:          {len(app_state_machine())}")
    print()
    print("  Guvenlik prensibi: Mobil uygulama tek basina guvenlik")
    print("  zincirinin parcasi DEGILDIR. Her kritik komut (ARM, START)")
    print("  donanim tarafinda ikinci faktor gerektirir.")
    print()

    # --- Figure ---
    fig_path = Path("zt008_5_screens.png")
    fig_ok = save_screens_figure(str(fig_path))
    print(f"Figure saved: {fig_path} ({'OK' if fig_ok else 'FAIL'})")
    print()

    out = Path("zt008_5_results.json")
    out.write_text(json.dumps({
        "stage": "ZT-008.5",
        "mode": a.mode,
        "tech_stack": tech_stack_options(),
        "screens": screen_list(),
        "ble_protocol": ble_protocol(),
        "wifi_protocol": wifi_protocol(),
        "app_states": list(app_state_machine().keys()),
    }, indent=2), encoding="utf-8")
    print(f"Results written to {out}")
    print("SCIENTIFIC STATUS: conceptual mobile spec only.")


if __name__ == "__main__":
    main()
'''

TST_CODE = r'''
"""ZT-008.5 tests."""
import numpy as np
from zt005.run_zt008_5 import (
    tech_stack_options, screen_list, ble_protocol,
    wifi_protocol, app_state_machine, compute_link_requirements,
)


def test_tech_stack_options():
    opts = tech_stack_options()
    assert len(opts) >= 4
    # Flutter should be PRIMARY
    flutter = [o for o in opts if "Flutter" in o["name"]]
    assert len(flutter) == 1
    assert flutter[0]["recommendation"] == "PRIMARY"


def test_screens_include_required():
    screens = [s["name"] for s in screen_list()]
    for required in ["Splash", "Home / Dashboard", "Run Monitor",
                     "Arm / Ready"]:
        assert required in screens


def test_screens_have_actions():
    for s in screen_list():
        assert "actions" in s
        assert len(s["actions"]) > 0


def test_ble_protocol_has_service_and_chars():
    ble = ble_protocol()
    assert "service_uuid" in ble
    assert "characteristics" in ble
    assert len(ble["characteristics"]) >= 3


def test_ble_command_char_lists_commands():
    ble = ble_protocol()
    cmd_char = [c for c in ble["characteristics"] if c["name"] == "Command"]
    assert len(cmd_char) == 1
    assert "ARM" in cmd_char[0]["commands"]
    assert "ABORT" in cmd_char[0]["commands"]


def test_wifi_has_rest_and_ws():
    wifi = wifi_protocol()
    assert "rest_endpoints" in wifi
    assert "websocket" in wifi
    assert len(wifi["rest_endpoints"]) >= 5


def test_app_state_machine_starts_disconnected():
    sm = app_state_machine()
    assert "DISCONNECTED" in sm


def test_app_state_transitions_valid():
    sm = app_state_machine()
    all_states = set(sm.keys())
    for state, trans in sm.items():
        for action, target in trans.items():
            assert target in all_states, \
                f"{state}->{action}->{target} invalid"


def test_run_monitor_can_use_usb():
    req = compute_link_requirements()
    assert "Run_Monitor" in req
    assert req["Run_Monitor"]["usb"] is True


def test_splash_only_uses_ble():
    req = compute_link_requirements()
    assert req["Splash"]["ble"] is True
    assert req["Splash"]["wifi"] is False
'''

if __name__ == "__main__":
    SRC.write_text(SRC_CODE, encoding="utf-8")
    print(f"[OK] Wrote {SRC}")
    TST.write_text(TST_CODE, encoding="utf-8")
    print(f"[OK] Wrote {TST}")
    print()
    print("Now run:")
    print("  python -m pytest tests\\test_zt008_5.py -v")
    print("  python -m zt005.run_zt008_5 --mode quick")