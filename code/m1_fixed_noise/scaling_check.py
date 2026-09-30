"""Independent integer replay of the asymptotic majorants."""

def verify(row):
    from m1_fermionic_source.check import exact, keys, need
    keys(row, 'parent_locations_q_degree accounting_q_degree local_fault_q_degree '
              'synthesis_accuracy_q_degree universal_synthesis_and_bookkeeping_log_degree '
              'levels synthesis_state_q_degree synthesis_accounting_q_degree selected_level '
              'full_library_constant error_output auxiliary_records local_failure')
    exact([row[k] for k in ('parent_locations_q_degree', 'accounting_q_degree', 'local_fault_q_degree',
                            'synthesis_accuracy_q_degree', 'universal_synthesis_and_bookkeeping_log_degree')],
          [4, 4, -1, -16, 6])
    need(type(row['levels']) is list and len(row['levels']) == 4, 'complete comparison of encoding levels')
    power = 1
    for level, bound in enumerate(row['levels'], 1):
        power *= 2
        exact(bound, dict(level=level, fault_power=power, state_q_degree=4-power, accounting_q_degree=8-power))
    exact(row['selected_level'], 4)
    exact(row['synthesis_state_q_degree'], 4-16)
    exact(row['synthesis_accounting_q_degree'], 4+4-16)
    need(row['levels'][-1]['accounting_q_degree'] < 0 and row['levels'][-2]['accounting_q_degree'] == 0,
         'four-level sufficient bound; no three-level convergence asserted')
    exact(row['full_library_constant'], 'symbolic; not the recovery pair count')
    exact(row['error_output'], 'decoded logical state, named records and abort')
    exact(row['auxiliary_records'], 'retained; no all-record closeness claim')
    exact(row['local_failure'], 'complete faulty output; only outermost failure aborts application')
