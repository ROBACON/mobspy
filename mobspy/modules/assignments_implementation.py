"""Implement SBML assignment rules that bind species counts to expressions."""

from __future__ import annotations

from contextvars import ContextVar
from re import compile as re_compile
from re import escape as re_escape
from typing import TYPE_CHECKING, Any

from mobspy.constants import ALL_CHAR, ASSIGNMENT_PREFIX, DOT_SEPARATOR
from mobspy.exceptions import CompilationError

if TYPE_CHECKING:
    from mobspy.types import AssignmentsForSbml

from mobspy.modules.mobspy_expressions import MobsPyExpression as mbe_MobsPyExpression
from mobspy.modules.species_string_generator import (
    construct_all_combinations as ssg_construct_all_combinations,
)
from mobspy.modules.species_string_generator import (
    construct_species_char_list as ssg_construct_species_char_list,
)
from mobspy.types import AssignmentData

_asg_context_cv: ContextVar[bool] = ContextVar("_asg_context_cv", default=False)


class Assignment_Operator:
    """Manages assignment context and compiles assignment expressions for SBML.

    Acts as a context manager that activates/deactivates assignment mode,
    and provides static methods for arithmetic on assignment expressions.
    """

    regex_pattern: str = r"\(\$arg(?:\.[^\s().]+)?\)(?=[^\s()]|$)"

    @staticmethod
    def find_arg_strings(input_string: str) -> list[str]:
        """Extract ``$``-prefixed argument tokens from an assignment expression."""
        arg_strings: list[str] = []
        flag_found = False
        stack = ""
        for char in input_string:
            if char == "$":
                flag_found = True

            if not flag_found:
                continue
            if char != ")":
                stack += char
            else:
                flag_found = False
                arg_strings.append(stack)
                stack = ""

        return arg_strings

    def __enter__(self) -> Assignment_Operator:
        _asg_context_cv.set(True)
        return self

    def set_context(self) -> None:
        """Activate the assignment context."""
        _asg_context_cv.set(True)

    def __exit__(self, *args: Any) -> None:
        _asg_context_cv.set(False)

    def reset_context(self) -> None:
        """Deactivate the assignment context."""
        _asg_context_cv.set(False)

    def check_context(self) -> bool:
        """Return whether the assignment context is currently active."""
        return _asg_context_cv.get()

    @staticmethod
    def check_arguments(
        first: Any,
        second: Any,
    ) -> tuple[mbe_MobsPyExpression, mbe_MobsPyExpression]:
        """Coerce both operands into MobsPyExpression instances.

        Prepares them for assignment arithmetic.
        """
        spe_list_first: list[Any] = []
        spe_list_second: list[Any] = []

        if hasattr(first, "get_spe_object"):
            spe_list_first = [first]

        if hasattr(second, "get_spe_object"):
            spe_list_second = [second]

        if isinstance(first, (int, float)):
            first = mbe_MobsPyExpression(
                str(first),
                species_object=None,
                dimension=None,
                count_in_model=True,
                concentration_in_model=False,
                count_in_expression=False,
                concentration_in_expression=False,
                species_list_operation_order=spe_list_first,
            )

        if isinstance(second, (int, float)):
            second = mbe_MobsPyExpression(
                str(second),
                species_object=None,
                dimension=None,
                count_in_model=True,
                concentration_in_model=False,
                count_in_expression=False,
                concentration_in_expression=False,
                species_list_operation_order=spe_list_second,
            )

        if not isinstance(first, mbe_MobsPyExpression):
            first = mbe_MobsPyExpression(
                "(" + ASSIGNMENT_PREFIX + str(first) + ")",
                species_object=None,
                dimension=None,
                count_in_model=True,
                concentration_in_model=False,
                count_in_expression=False,
                concentration_in_expression=False,
                species_list_operation_order=spe_list_first,
            )
        if not isinstance(second, mbe_MobsPyExpression):
            second = mbe_MobsPyExpression(
                "($asg_" + str(second) + ")",
                species_object=None,
                dimension=None,
                count_in_model=True,
                concentration_in_model=False,
                count_in_expression=False,
                concentration_in_expression=False,
                species_list_operation_order=spe_list_second,
            )
        return first, second

    @staticmethod
    def add(first: Any, second: Any) -> mbe_MobsPyExpression:
        """Build an addition expression from two assignment operands."""
        first, second = Assignment_Operator.check_arguments(first, second)
        return first + second  # type: ignore[no-any-return]

    @staticmethod
    def sub(first: Any, second: Any) -> mbe_MobsPyExpression:
        """Build a subtraction expression from two assignment operands."""
        first, second = Assignment_Operator.check_arguments(first, second)
        return first - second  # type: ignore[no-any-return]

    @staticmethod
    def mul(first: Any, second: Any) -> mbe_MobsPyExpression:
        """Build a multiplication expression from two assignment operands."""
        first, second = Assignment_Operator.check_arguments(first, second)
        return first * second  # type: ignore[no-any-return]

    @staticmethod
    def div(first: Any, second: Any) -> mbe_MobsPyExpression:
        """Build a division expression from two assignment operands."""
        first, second = Assignment_Operator.check_arguments(first, second)
        return first / second  # type: ignore[no-any-return]

    @staticmethod
    def pow(first: Any, second: Any) -> mbe_MobsPyExpression:
        """Build an exponentiation expression from two assignment operands."""
        first, second = Assignment_Operator.check_arguments(first, second)
        return first**second  # type: ignore[no-any-return]

    @staticmethod
    def generate_replacement_in_expression(
        express_spe: str,
        ortogonal_vector_structure: dict[str, Any],
        meta_species_in_model: list[Any],
        expression_tuple: tuple[Any, str],
    ) -> str:
        """Expand a meta-species token into a sum of concrete species.

        Returns a string of all matching concrete species joined by '+'.
        """
        spe_str_raw = express_spe.replace(ASSIGNMENT_PREFIX, "")
        spe_str = spe_str_raw.split(".")

        # CHECK HERE FOR MISSING SPECIES IN MODEL
        for meta_spe in meta_species_in_model:
            if spe_str[0] == str(meta_spe):
                spe_object = meta_spe
                break
        else:
            spe_name = str(expression_tuple[0][0])
            error_message = (
                f"Assignment {spe_name}, {expression_tuple[0][1]}: "
                f"{expression_tuple[1].replace(ASSIGNMENT_PREFIX, '')} failed\n"
                "One of the meta-species in the assignment"
                " expression was not found in the model"
            )
            raise CompilationError(error_message)

        if len(spe_str) == 1:
            str_comb = ssg_construct_all_combinations(
                spe_object, set(), ortogonal_vector_structure, DOT_SEPARATOR
            )
        else:
            str_comb = ssg_construct_all_combinations(
                spe_object, set(spe_str[1:]), ortogonal_vector_structure, DOT_SEPARATOR
            )
        str_comb.sort()

        to_replace: str = ""
        for i, e in enumerate(str_comb):
            if i == 0:
                to_replace = e
            else:
                to_replace += "+" + e
        return to_replace

    @staticmethod
    def process_assignments(
        asg_expression: str,
        ortogonal_vector_structure: dict[str, Any],
        meta_species_in_model: list[Any],
        for_error_tuple: tuple[Any, str],
    ) -> str:
        """Replace all meta-species tokens in an assignment expression.

        Substitutes concrete species sums for each token.
        """
        spe_in_expression = Assignment_Operator.find_arg_strings(asg_expression)
        replacement_dict: dict[str, str] = {}
        for spe in spe_in_expression:
            replacement_dict[spe] = (
                Assignment_Operator.generate_replacement_in_expression(
                    spe,
                    ortogonal_vector_structure,
                    meta_species_in_model,
                    for_error_tuple,
                )
            )

        for key, item in replacement_dict.items():
            asg_expression = Assignment_Operator.replace_agn_expr(
                asg_expression, key, item
            )

        return asg_expression

    @staticmethod
    def compile_assignments_for_sbml(
        unprocessed_asgns: dict[Any, Any],
        ortogonal_vector_structure: dict[str, Any],
        meta_species_in_model: Any,
    ) -> AssignmentsForSbml:
        """Compile raw assignment definitions into SBML-ready AssignmentData entries."""
        assignments_for_sbml: AssignmentsForSbml = {}
        assignment_counter = 0
        for asg in unprocessed_asgns:
            if ALL_CHAR not in asg[1]:
                continue

            spe_to_asgn = ssg_construct_all_combinations(
                asg[0], asg[1], ortogonal_vector_structure, symbol=DOT_SEPARATOR
            )

            asgn_expression = Assignment_Operator.process_assignments(
                str(unprocessed_asgns[asg]),
                ortogonal_vector_structure,
                meta_species_in_model,
                (asg, str(unprocessed_asgns[asg])),
            )

            for spe in spe_to_asgn:
                key = "assignment_" + str(assignment_counter)
                assignments_for_sbml[key] = AssignmentData(
                    species=spe,
                    expression=asgn_expression,
                )
                assignment_counter += 1

        for asg in unprocessed_asgns:
            if ALL_CHAR in asg[1]:
                continue

            spe_to_asgn_result = ssg_construct_species_char_list(
                asg[0], asg[1], ortogonal_vector_structure, symbol=DOT_SEPARATOR
            )

            asgn_expression = Assignment_Operator.process_assignments(
                str(unprocessed_asgns[asg]),
                ortogonal_vector_structure,
                meta_species_in_model,
                (asg, str(unprocessed_asgns[asg])),
            )

            key = "assignment_" + str(assignment_counter)
            assignments_for_sbml[key] = AssignmentData(
                species=str(spe_to_asgn_result),
                expression=asgn_expression,
            )
            assignment_counter += 1

        return assignments_for_sbml

    @staticmethod
    def replace_agn_expr(assignment_ex: str, to_replace: str, replacement: str) -> str:
        """Replace a token in an assignment expression.

        Uses word-boundary-aware regex to avoid partial matches.
        """
        pattern = re_compile(re_escape(to_replace) + r"(?![a-zA-Z0-9_])")
        return pattern.sub(replacement, assignment_ex)


Assign = Assignment_Operator()


class Asg:
    """Captures an assignment target for deferred evaluation.

    Stores species and characteristics.
    """

    assignments: dict[Any, Any] = {}

    def __init__(self, meta_spe: Any, species_or_reacting: bool) -> None:
        Assign.set_context()
        self.meta_spe: list[Any] = []
        self.asgn_key: list[tuple[Any, tuple[Any, ...]]] = []
        if species_or_reacting:
            self.meta_spe.append(meta_spe)
            self.asgn_key.append((meta_spe, ()))
        else:
            for reacting_spe in meta_spe.list_of_reactants:
                self.meta_spe.append(reacting_spe["object"])
                self.asgn_key.append(
                    (reacting_spe["object"], tuple(reacting_spe["characteristics"]))
                )
        self.species_or_reacting = species_or_reacting

    def __call__(self, assignment: Any) -> None:
        for spe, key in zip(self.meta_spe, self.asgn_key):
            spe._assignments[key] = assignment
        Assign.reset_context()

    def __getattr__(self, item: str) -> None:
        raise CompilationError(
            "Assignments must be the last query in the"
            " stack - Ex: A.young.blue.assign()"
        )
