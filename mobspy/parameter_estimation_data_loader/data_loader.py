"""Load experimental time-series data into structures suitable for parameter fitting."""

from __future__ import annotations

from typing import Any

import mobspy.data_handler.time_series_object as tso
from mobspy.exceptions import ValidationError


class Experimental_Data_Holder:
    """Mixin that stores experimental data for parameter estimation."""

    def __init__(self) -> None:
        self.experimental_data: Any = None

    def load_experiment_data(self, data: Any) -> None:
        """Store experimental data for parameter estimation.

        Args:
            data: List of dicts or a MobsPyList_of_TS result.

        Raises:
            ValidationError: If the data format is invalid.
        """
        flag_jump_checks = True if isinstance(data, tso.MobsPyList_of_TS) else False  # noqa: SIM210

        if type(data) != list and not flag_jump_checks:  # noqa: E721
            raise ValidationError(
                "Data added must be in the format of list with"
                " each element being a dictionary "
                "with species names and time as keys"
                " or a MobsPy results object"
            )
        for e in data:
            if type(e) != dict and not flag_jump_checks:  # noqa: E721
                raise ValidationError(
                    "Data added must be in the format of list"
                    " with each element being a dictionary "
                    "with species names and time as keys"
                    " or a MobsPy results object"
                )

        self.experimental_data = data
