"""
The any_species.py model is responsible
for defining the Context_specie_named_any class.
"""

from __future__ import annotations

from typing import Any as TypingAny
from typing import NoReturn

from mobspy.constants import CONTEXT_ANY_SPECIES_NAME
from mobspy.exceptions import ValidationError
from mobspy.modules.species import Species


def _get_any_chars() -> set[str]:
    """Return the current Any characteristics set from session context."""
    from mobspy.modules.session_context import get_session  # noqa: PLC0415

    return get_session().any_chars


def _get_any_stack() -> list[set[str]]:
    """Return the nested Any context stack from session context."""
    from mobspy.modules.session_context import get_session  # noqa: PLC0415

    return get_session().any_stack


class Context_specie_named_any(Species):  # noqa: N801
    """
    Class which inherits from Species. It only has one
    object, Any, which is defined at the end of this
    script. It is used to simplify the syntax of reactions
    and count setting inside the body of a
    "with Any.example_characteristic1 :" statement.
    Is is compatible with event_condition and event_time
    methods of the Simulation class.
    It can be nested and used several time in the same
    with statement.
    """

    def __getattr__(self, item: str) -> TypingAny:  # type: ignore[override]
        """
        This method is called when an attribute is called
        on the Any specie. It is used to add the
        characteristic to the currently active Any context.

        Args:
            item: Characteristic to be added to the currently active Any context.
        """
        if item.startswith("_"):
            raise AttributeError(item)

        from mobspy.modules.session_context import get_session  # noqa: PLC0415

        get_session().any_building = True
        _get_any_chars().add(item)
        return self

    def __enter__(self) -> int:
        """
        Context manager for Any's characteristics.
        Called in "with Any.example_characteristic :" format, when entering.
        """
        from mobspy.modules.session_context import get_session  # noqa: PLC0415

        session = get_session()
        if not session.any_building:
            raise ValidationError(
                "Characteristics cannot be added to the Any specie outside of a context"
            )
        session.any_building = False
        self.context_initiator_for_meta_specie_named_any()
        return 0

    def __exit__(self, *args: TypingAny) -> None:
        """
        Context manager for Any's characteristics.
        Called in "with Any.example_characteristic :" format, when exiting.
        """
        self.context_finish_for_meta_specie_named_any()

    def context_initiator_for_meta_specie_named_any(self) -> None:
        """
        This adds the current context in the nested stack,
        and then updates the Any context in all meta-species.
        """
        chars = _get_any_chars()
        _get_any_stack().append(set(chars))
        Species.update_meta_specie_named_any_context(
            Species.get_meta_specie_named_any_context().union(chars)
        )

    def context_finish_for_meta_specie_named_any(self) -> None:
        """
        This removes the context which is ending from
        the nested stack and updates the current Any context.
        Then, it updates the Any context in all meta-species.
        """
        from mobspy.modules.session_context import get_session  # noqa: PLC0415

        stack = _get_any_stack()
        previous_chars = stack.pop()
        session = get_session()
        if len(stack) > 0:
            session.any_chars = stack[-1]
        else:
            session.any_chars = set()
        Species.update_meta_specie_named_any_context(
            Species.get_meta_specie_named_any_context() - previous_chars
        )

    def __call__(self, quantity: TypingAny) -> None:  # noqa: ARG002
        """
        The call operator is overloaded as the Any specie
        cannot be called, as it is not supposed to be used
        this way. Thus, it raises an error when called.
        """
        raise ValidationError("The Any species cannot be called")

    def __add__(self, other: TypingAny) -> None:
        """
        The add operator is overloaded as the Any species
        cannot be added. Raises an error when added.
        """
        raise ValidationError("The Any species cannot be added")

    def __radd__(self, other: TypingAny) -> None:
        """
        The add operator is overloaded as the Any species
        cannot be added. Raises an error when added.
        """
        raise ValidationError("The Any species cannot be added")

    def __rmul__(self, other: TypingAny) -> None:
        """
        The multiplication operator is overloaded as the
        Any species cannot be multiplied. Raises an error
        when multiplied.
        """
        raise ValidationError("The Any species cannot be multiplied")

    def __rshift__(self, other: TypingAny) -> NoReturn:
        """
        The >> operator is overloaded so that  it raises an error when used.
        """
        raise ValidationError("The >> operator cannot be used on the Any species")


# Any is the only object of the Any_specie class that will be used. It is defined here.
__SAny = Context_specie_named_any(CONTEXT_ANY_SPECIES_NAME)
Any = __SAny
