"""Replay the actual finite trial and endogenous history evidence."""
import copy
import hashlib
import json

import pytest

import source_derived_public_checkpoint_packet as source


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def rehash(packet, manifest):
    manifest.pop("terminal_fiber_manifest_sha256", None)
    manifest["terminal_fiber_manifest_sha256"] = digest(manifest)
    packet["terminal_fiber_manifest_hash"] = manifest["terminal_fiber_manifest_sha256"]
    packet.pop("packet_sha256", None)
    packet["packet_sha256"] = digest(packet)


@pytest.fixture
def evidence():
    manifest = source.build_terminal_fiber_manifest()
    return source.build_source_derived_packet(manifest), manifest


def test_public_rejects_empty_actual_terminal_fiber(evidence):
    packet, manifest = evidence
    manifest.update(trials=[], trial_count=0, terminal_world_ids=[])
    rehash(packet, manifest)
    assert source.certify_source_derived_packet(packet, terminal_manifest=manifest)["status"] != "PASS"


def test_public_rejects_actual_target_reading_history(evidence):
    packet, manifest = evidence
    for sid, history in packet["semantic_histories"].items():
        history.update(events=["read:expected_answer:24"], semantic_event_count=1,
                       uses_external_target=True)
        packet["reachability_witnesses"][sid] = history["events"]
    rehash(packet, manifest)
    assert source.certify_source_derived_packet(packet, terminal_manifest=manifest)["status"] != "PASS"


@pytest.mark.parametrize("fault", ["missing_trial", "duplicate_trial", "wrong_candidate",
                                  "wrong_membership", "wrong_terminal_ids", "nonterminal_packet",
                                  "bad_constructor", "bad_self_hash", "incomplete_claim"])
def test_manifest_replays_actual_trials_and_membership(evidence, fault):
    packet, manifest = evidence
    if fault == "missing_trial":
        manifest["trials"].pop()
    elif fault == "duplicate_trial":
        manifest["trials"][-1] = copy.deepcopy(manifest["trials"][0])
    elif fault == "wrong_candidate":
        row = manifest["trials"][1]
        row["candidate"]["edges"] = copy.deepcopy(manifest["trials"][0]["candidate"]["edges"])
        row["candidate_sha256"] = digest(row["candidate"])
    elif fault == "wrong_membership":
        manifest["trials"][1]["terminal"] = True
    elif fault == "wrong_terminal_ids":
        manifest["terminal_world_ids"] = []
    elif fault == "nonterminal_packet":
        packet["terminal_id"] = manifest["trials"][1]["trial_id"]
    elif fault == "bad_constructor":
        manifest["trial_universe_constructor_sha256"] = "f" * 64
        packet["trial_universe_constructor_hash"] = "f" * 64
    elif fault == "incomplete_claim":
        manifest["terminal_fiber_completeness_certificate"]["complete"] = False
    rehash(packet, manifest)
    if fault == "bad_self_hash":
        manifest["terminal_fiber_manifest_sha256"] = "f" * 64
        packet["terminal_fiber_manifest_hash"] = "f" * 64
        packet.pop("packet_sha256")
        packet["packet_sha256"] = digest(packet)
    assert source._verify_terminal_manifest(packet, manifest)["status"] != "PASS"
    assert source.certify_source_derived_packet(packet, terminal_manifest=manifest)["status"] != "PASS"


def test_reordered_trial_and_set_presentations_remain_valid(evidence):
    packet, manifest = evidence
    manifest["trials"].reverse()
    for row in manifest["trials"]:
        candidate = row["candidate"]
        candidate["edges"] = [list(reversed(edge)) for edge in reversed(candidate["edges"])]
        candidate["oriented_slots"].reverse()
        row["candidate_sha256"] = digest(candidate)
    rehash(packet, manifest)
    assert source._verify_terminal_manifest(packet, manifest)["status"] == "PASS"
    receipt = source.certify_source_derived_packet(packet, terminal_manifest=manifest)
    assert receipt["status"] == "PASS"
    assert receipt["capacity"]["robust_closure_at_frozen_D"] is True


@pytest.mark.parametrize("fault", ["omitted_trial", "wrong_materialization"])
def test_shared_constructor_bug_cannot_define_its_own_oracle(monkeypatch, fault):
    if fault == "omitted_trial":
        original = source._trial_specs
        def omitted_trial():
            return original()[:-1]
        monkeypatch.setattr(source, "_trial_specs", omitted_trial)
    else:
        original = source._materialize_trial
        def wrong_materialization(spec):
            candidate = original(spec)
            if spec["mutation_kind"] == "delete_edge":
                # Keep the right number of faults and terminal/nonterminal
                # labels, but remove a different edge than the named input.
                chosen = source.icosahedral_edges()[-1]
                candidate["edges"] = [list(edge) for edge in source.icosahedral_edges()
                                      if edge != chosen]
            return candidate
        monkeypatch.setattr(source, "_materialize_trial", wrong_materialization)
    manifest = source.build_terminal_fiber_manifest()
    if fault == "omitted_trial":
        # Make every claimed count internally consistent with the incomplete
        # enumeration, so only independent knowledge of all 67 faults detects it.
        manifest["terminal_fiber_completeness_certificate"]["expected_trial_count"] = manifest["trial_count"]
    packet = source.build_source_derived_packet(manifest)
    rehash(packet, manifest)
    assert source._verify_terminal_manifest(packet, manifest)["status"] != "PASS"
    assert source.certify_source_derived_packet(packet, terminal_manifest=manifest)["status"] != "PASS"


@pytest.mark.parametrize("fault", ["empty", "missing", "extra", "reversed", "target_flag",
                                  "executor_flag", "wrong_count", "wrong_source", "disagreeing_witness",
                                  "unreached_parent", "nonedge", "duplicate_child", "wrong_slot"])
def test_history_replay_rejects_missing_or_impossible_custody(evidence, fault):
    packet, manifest = evidence
    histories = packet["semantic_histories"]
    sid = packet["public_section_aliases"]["north/write"]
    history = histories[sid]
    history["events"] = list(history["events"])
    events = history["events"]
    if fault == "empty":
        histories.clear()
    elif fault == "missing":
        histories.pop(sid)
    elif fault == "extra":
        histories["unbound_record"] = copy.deepcopy(history)
    elif fault == "reversed":
        events.reverse()
    elif fault == "target_flag":
        history["uses_external_target"] = True
    elif fault == "executor_flag":
        history["uses_executor_metadata"] = 0
    elif fault == "wrong_count":
        history["semantic_event_count"] += 1
    elif fault == "wrong_source":
        history["source_slot"] = "south/write"
    elif fault == "disagreeing_witness":
        packet["reachability_witnesses"][sid] = list(reversed(events))
    elif fault == "unreached_parent":
        events[2] = "repair-propagate:upper_0->lower_0:north/write"
    elif fault == "nonedge":
        events[2] = "repair-propagate:north->south:north/write"
    elif fault == "duplicate_child":
        events[3] = events[2]
    elif fault == "wrong_slot":
        events[2] = events[2].replace("north/write", "north/check")
    if fault != "disagreeing_witness":
        packet["reachability_witnesses"][sid] = list(events)
    rehash(packet, manifest)
    assert source._verify_source_histories(packet)["status"] != "PASS"
    assert source.certify_source_derived_packet(packet, terminal_manifest=manifest)["status"] != "PASS"


def test_alternative_valid_spanning_tree_execution_remains_valid(evidence):
    packet, manifest = evidence
    adjacency = {port: set() for port in source.PORTS}
    for left, right in source.icosahedral_edges():
        adjacency[left].add(right)
        adjacency[right].add(left)
    for sid, history in packet["semantic_histories"].items():
        slot, root = history["source_slot"], history["source_port"]
        reached, propagated = {root}, []
        def visit(parent):
            for child in sorted(adjacency[parent], reverse=True):
                if child not in reached:
                    reached.add(child)
                    propagated.append(f"repair-propagate:{parent}->{child}:{slot}")
                    visit(child)
        visit(root)
        assert len(reached) == 12 and len(propagated) == 11
        events = history["events"][:2] + propagated + history["events"][-2:]
        assert events != history["events"]
        history["events"] = events
        packet["reachability_witnesses"][sid] = list(events)
    rehash(packet, manifest)
    assert source._verify_source_histories(packet)["status"] == "PASS"
    assert source.certify_source_derived_packet(packet, terminal_manifest=manifest)["status"] == "PASS"


@pytest.mark.parametrize("fault", ["empty_fiber", "target_history"])
def test_serialized_source_refusal_cannot_feed_direct_n_closure(evidence, fault, tmp_path, monkeypatch):
    import direct_n_closure_verdict as downstream

    packet, manifest = evidence
    if fault == "empty_fiber":
        manifest.update(trials=[], trial_count=0, terminal_world_ids=[])
    else:
        sid = next(iter(packet["semantic_histories"]))
        packet["semantic_histories"][sid]["uses_external_target"] = True
    rehash(packet, manifest)
    receipt = source.certify_source_derived_packet(packet, terminal_manifest=manifest)
    path = tmp_path / "fixed-source-certificate.json"
    path.write_text(json.dumps(receipt, allow_nan=False), encoding="utf-8")
    monkeypatch.setattr(downstream, "FIXED_CERTIFICATE_PATH", path)
    with pytest.raises(ValueError, match="fixed-cutoff parent is not attained"):
        downstream.build_verdict()
