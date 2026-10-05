
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
