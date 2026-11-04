"""The recipe registry: the single description of every recipe.

Each recipe module defines its recipes in a ``*_RECIPES_SCHEMAS``
dictionary (its brief and details, its parameters with their types,
defaults and descriptions, and its flags). The registry reads them once,
at import, into frozen :class:`RecipeSpec` and :class:`ParamSpec` objects;
the GUI (menus, dialogs, requirements, branching), the scripts and the MCP
server read the registry, never the dictionaries, so that a recipe is
described in one place and nothing changes its description at run time
(findings register F12 and F79; ``ROADMAP.md``, phase 6-bis).

Adding a recipe means adding its entry to the dictionary of its module
(and its method to the recipe class): ``tests/test_recipe_registry.py``
checks that the two agree.

Keys of a recipe entry:

``brief``
    One line, the name of the recipe in the menus and dialogs.
``details``
    A longer description (optional).
``url``
    The page of the documentation (optional).
``params``
    The list of the parameters (keys below).
``gui_hidden``
    The recipe is not in the menus (it is run by the GUI or by scripts).
``branching``
    The recipe opens a new session (or several) instead of giving the next
    state of the same session.
``requires``
    What the recipe needs from the session before it runs: ``"z_em"``,
    ``"resolution"``. Today the GUI asks for them; phase 6-bis extends
    this to the other prerequisites ("Prerequisites of the recipes").

Keys of a parameter: ``name``, ``type``, ``default``, ``doc`` (optional,
the name if missing), ``gui_hidden`` (not shown in the dialog), and for the
conversion of its value (:mod:`astrocook.recipes.params`): ``optional``
(None is a value: no value, or automatic), ``keep`` (special strings passed
as they are, e.g. ``'_current_'``), ``min`` and ``max`` (bounds of a
number, included), ``gt`` (a bound excluded: the value must be above it),
``allow_nan`` (NaN is a value; infinities never are).
"""
from __future__ import annotations

import copy
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any, Dict, Iterator, Mapping, Optional, Tuple

#: The categories, in the order of the menus
CATEGORIES: Tuple[str, ...] = ('file', 'edit', 'flux', 'continuum', 'absorbers')

#: What a recipe can require (the key ``requires``)
REQUIREMENTS = frozenset({'z_em', 'resolution'})

_RECIPE_KEYS = {'brief', 'details', 'url', 'params', 'gui_hidden', 'branching',
                'requires'}
_PARAM_KEYS = {'name', 'type', 'default', 'doc', 'gui_hidden', 'optional', 'keep',
               'min', 'max', 'gt', 'allow_nan'}


@dataclass(frozen=True)
class ParamSpec:
    """A parameter of a recipe."""
    #: Name, as in the signature of the method
    name: str
    #: Python type of the value (``str``, ``float``, ``int``, ``bool``, ``list``)
    type: type
    #: Description, shown by the dialogs and the MCP tools
    doc: str
    #: Not shown in the dialog of the recipe
    gui_hidden: bool = False
    #: None is a value (no value, or automatic)
    optional: bool = False
    #: Special strings passed as they are (e.g. ``'_current_'``)
    keep: Tuple[str, ...] = ()
    #: Bounds of a number (None: no bound): ``min`` and ``max`` included,
    #: ``gt`` excluded (the value must be above it)
    min: Optional[float] = None
    max: Optional[float] = None
    gt: Optional[float] = None
    #: NaN is a value of a number
    allow_nan: bool = False
    _default: Any = field(default=None, repr=False)

    @property
    def default(self) -> Any:
        """The default value (a copy: a list default cannot be changed)."""
        return copy.deepcopy(self._default)


@dataclass(frozen=True)
class RecipeSpec:
    """A recipe: what it is, its parameters and its flags."""
    #: Menu of the recipe (``file``, ``edit``, ``flux``, ``continuum``, ``absorbers``)
    category: str
    #: Name of the method of the recipe class
    name: str
    #: One line, the name of the recipe in the menus and dialogs
    brief: str
    #: A longer description
    details: str = ""
    #: The page of the documentation
    url: Optional[str] = None
    #: The parameters, in the order of the dialog
    params: Tuple[ParamSpec, ...] = ()
    #: Not in the menus
    gui_hidden: bool = False
    #: Opens a new session instead of giving the next state of the same one
    branching: bool = False
    #: What the recipe needs from the session (a subset of :data:`REQUIREMENTS`)
    requires: frozenset = frozenset()

    def param(self, name: str) -> ParamSpec:
        """The parameter called ``name`` (``KeyError`` if there is none)."""
        for p in self.params:
            if p.name == name:
                return p
        raise KeyError(f"{self.category}.{self.name} has no parameter {name!r}")

    def defaults(self) -> Dict[str, Any]:
        """The defaults of all the parameters, by name."""
        return {p.name: p.default for p in self.params}


def _param_spec(where: str, entry: Mapping[str, Any]) -> ParamSpec:
    unknown = set(entry) - _PARAM_KEYS
    if unknown:
        raise ValueError(f"{where}: unknown keys {sorted(unknown)}")
    return ParamSpec(name=entry['name'], type=entry['type'],
                     doc=entry.get('doc', entry['name']),
                     gui_hidden=bool(entry.get('gui_hidden', False)),
                     optional=bool(entry.get('optional', False)),
                     keep=tuple(entry.get('keep', ())),
                     min=entry.get('min'), max=entry.get('max'), gt=entry.get('gt'),
                     allow_nan=bool(entry.get('allow_nan', False)),
                     _default=copy.deepcopy(entry['default']))


def build_spec(category: str, name: str, entry: Mapping[str, Any]) -> RecipeSpec:
    """The :class:`RecipeSpec` of one entry of a ``*_RECIPES_SCHEMAS``
    dictionary; ``ValueError`` on an unknown key or requirement."""
    where = f"{category}.{name}"
    unknown = set(entry) - _RECIPE_KEYS
    if unknown:
        raise ValueError(f"{where}: unknown keys {sorted(unknown)}")
    requires = frozenset(entry.get('requires', ()))
    if requires - REQUIREMENTS:
        raise ValueError(f"{where}: unknown requirements {sorted(requires - REQUIREMENTS)}")
    params = tuple(_param_spec(f"{where}", p) for p in entry.get('params', []))
    return RecipeSpec(category=category, name=name, brief=entry['brief'],
                      details=entry.get('details', ""), url=entry.get('url'),
                      params=params, gui_hidden=bool(entry.get('gui_hidden', False)),
                      branching=bool(entry.get('branching', False)),
                      requires=requires)


def _build() -> Mapping[str, Mapping[str, RecipeSpec]]:
    from astrocook.recipes.absorbers import ABSORBERS_RECIPES_SCHEMAS
    from astrocook.recipes.continuum import CONTINUUM_RECIPES_SCHEMAS
    from astrocook.recipes.edit import EDIT_RECIPES_SCHEMAS
    from astrocook.recipes.file import FILE_RECIPES_SCHEMAS
    from astrocook.recipes.flux import FLUX_RECIPES_SCHEMAS
    sources = {'file': FILE_RECIPES_SCHEMAS, 'edit': EDIT_RECIPES_SCHEMAS,
               'flux': FLUX_RECIPES_SCHEMAS, 'continuum': CONTINUUM_RECIPES_SCHEMAS,
               'absorbers': ABSORBERS_RECIPES_SCHEMAS}
    registry = {}
    seen = {}
    for category in CATEGORIES:
        specs = {}
        for name, entry in sources[category].items():
            # The GUI and the log find a recipe by its name alone
            if name in seen:
                raise ValueError(f"Recipe {name!r} in both {seen[name]} and {category}")
            seen[name] = category
            specs[name] = build_spec(category, name, entry)
        registry[category] = MappingProxyType(specs)
    return MappingProxyType(registry)


#: The registry: category -> recipe name -> :class:`RecipeSpec` (read only)
RECIPES: Mapping[str, Mapping[str, RecipeSpec]] = _build()


def get_recipe(category: str, name: str) -> RecipeSpec:
    """The :class:`RecipeSpec` of a recipe (``KeyError`` if there is none)."""
    try:
        return RECIPES[category][name]
    except KeyError:
        raise KeyError(f"There is no recipe {category}.{name}") from None


def find_recipe(name: str) -> Optional[RecipeSpec]:
    """The :class:`RecipeSpec` of a recipe by its name alone (names are
    unique across the categories), or None."""
    for specs in RECIPES.values():
        if name in specs:
            return specs[name]
    return None


def all_recipes() -> Iterator[RecipeSpec]:
    """All the recipes, by category in the order of the menus."""
    for category in CATEGORIES:
        yield from RECIPES[category].values()
