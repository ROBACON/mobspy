from test_script import *

# This is here because pytest is too slow - so it's a faster test that avoid the collection step.
# Pytest is used in the GitHub while pushing
# This script is oriented towards fast debugging and development
# Also in the future, it might be nice to add an automated way to get the test names
test_list = [test_model_1, test_model_2, test_model_3, test_model_4, test_model_5, test_model_6, test_model_7,
             test_orthogonal_spaces, test_average_value, test_hybrid_sim, test_concatenated_simulation,
             test_event_type, test_reacting_species_event, test_unit_event_test, test_reaction_deactivation,
             test_double_rate, test_single_rate, test_triple_rate, test_stochastic_event_duration,
             test_logic_operator_syntax, test_stack_position, test_empty_arguments,
             test_conditional_between_meta_species, test_conditional_between_meta_species_2,
             test_event_reaction_not_allowed, all_test, all_test_2, test_error_mult, test_set_counts,
             test_bool_error, test_event_all, test_one_value_concatenation_sim, test_crash_after_modification,
             test_unit_bi_dimension, test_bi_dimensional_rates, test_dimension_in_function_only,
             test_multiple_simulation_counts, test_string_events_assignment, test_plotting,
             test_volume_after_sim, test_parameters_with_sbml, test_shared_parameter_name,
             test_set_counts_parameters, test_repeated_parameters, initial_expression_test,
             test_wrong_dimension_error, test_more_than_used, zero_rate_test, test_wrong_rate,
             test_conversion_outside, test_first_characteristic_in_reacting_species, test_model_reference,
             test_sbml_generation, test_multi_sim_sbml, test_inline_comment,
             test_with_statement_any_and_species_characteristics, test_with_statement_on_any_and_event,
             test_matching_characteristic_rate, test_changes_after_compilation, test_proper_unit_context_exit,
             test_run_args, test_unit_args, test_multi_parameters_in_run, test_output_concentration_in_multi_sim,
             test_parameter_operation_in_rate, test_multi_parameter_with_expression, test_double_parameters_with_units,
             test_parameters_with_units, test_convert_back_parameter,
             test_numpy_in_expression_function, test_numpy_in_rates,
             test_numpy_in_counts, test_numpy_in_set_counts, test_multi_methods_plot, test_unit_x_conversion,
             test_Silicon_valley, test_replacing_species_name_in_expression, test_basic_assignment,
             test_all_asgn_ops, text_complex_assignments,
             text_assign_context_exit, text_even_more_complex_assignments, test_assign_context_complex,
             test_assign_context_constant, test_duration_with_run, test_rev, test_dimensionless_count,
             test_assignment_similar_species, test_blocked_names, test_blocked_names_2,
             test_update_parameter_for_multi_model, test_update_parameter_through_str,
             test_update_multiple_parameters_in_expression, test_update_parameter_with_unit,
             test_species_value_modification, test_all_value_modification, test_new_reversible_reaction_notation,
             test_2D_reaction_with_units, test_parameters_as_initial_values, test_parameters_in_lambda_expression]

# test_no_species_in_asg
# test_illegal_unit_op_in_assignment
# temporary_test_removal = [test_parameter_fit_with_units, test_multiple_runs_fit, test_simple_fit]
test_remov = [test_antimony_compose_model_gen, test_antimony_model]
sub_test = test_list
# sub_test = [test_parameters_with_sbml]
def perform_tests():
    any_failed = False
    for test in sub_test:
        try:
            test()
            print(f'Test {test} passed')
        except:
            print('\033[91m' + f'Test {test} failed' + '\033[0m', file=sys.stderr)
            any_failed = True
    if any_failed:
        assert False
perform_tests()