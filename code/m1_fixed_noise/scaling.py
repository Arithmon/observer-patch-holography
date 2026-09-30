"""Integer exponents of the proved majorants; no invented gadget threshold."""


def evidence():
    return dict(parent_locations_q_degree=4, accounting_q_degree=4,
                local_fault_q_degree=-1, synthesis_accuracy_q_degree=-16,
                universal_synthesis_and_bookkeeping_log_degree=6,
                levels=[dict(level=k, fault_power=2**k, state_q_degree=4-2**k,
                             accounting_q_degree=8-2**k) for k in range(1, 5)],
                synthesis_state_q_degree=-12, synthesis_accounting_q_degree=-8,
                selected_level=4, full_library_constant='symbolic; not the recovery pair count',
                error_output='decoded logical state, named records and abort',
                auxiliary_records='retained; no all-record closeness claim',
                local_failure='complete faulty output; only outermost failure aborts application')
