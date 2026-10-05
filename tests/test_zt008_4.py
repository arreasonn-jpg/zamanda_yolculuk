
"""ZT-008.4 tests - sensor + telemetry."""
import numpy as np
from zt005.run_zt008_4 import (
    sensor_list, telemetry_budget, channel_bandwidth,
    logging_capacity, time_sync_analysis,
)


def test_sensor_list_nonempty():
    assert len(sensor_list()) >= 6


def test_all_sensors_have_required_fields():
    for s in sensor_list():
        for key in ["name", "quantity", "range", "sample_rate_Hz",
                    "quantity_bytes"]:
            assert key in s


def test_telemetry_budget_positive():
    b = telemetry_budget()
    assert b["raw_bytes_per_s"] > 0
    assert b["compressed_bytes_per_s"] > 0
    assert b["compressed_bytes_per_s"] < b["raw_bytes_per_s"]


def test_bandwidth_links_ordered():
    bw = channel_bandwidth()
    assert bw["BLE 5.0 (effective)"] < bw["WiFi 802.11n (effective)"]
    assert bw["WiFi 802.11n (effective)"] < bw["USB 3.0 (effective)"]


def test_logging_capacity_scales():
    c1 = logging_capacity(1.0, 1e6)
    c2 = logging_capacity(2.0, 1e6)
    assert abs(c2 / c1 - 2.0) < 1e-9


def test_time_sync_has_csac():
    ts = time_sync_analysis()
    assert "atomic_portable" in ts


def test_csac_accuracy_is_best():
    ts = time_sync_analysis()
    csac_drift = ts["atomic_portable"]["drift_per_hour_ns"]
    gnss_drift = ts["gnss_disciplined"]["drift_per_hour_ns"]
    assert csac_drift < gnss_drift


def test_wifi_can_carry_raw_data():
    b = telemetry_budget()
    bw = channel_bandwidth()
    assert bw["WiFi 802.11n (effective)"] > b["raw_bytes_per_s"]


def test_ble_cannot_carry_raw_data():
    b = telemetry_budget()
    bw = channel_bandwidth()
    assert bw["BLE 5.0 (effective)"] < b["raw_bytes_per_s"]
