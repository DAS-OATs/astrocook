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
    ``"resolution"``, ``"continuum"`` (:mod:`astrocook.recipes.prerequisites`
    says what is missing and fills it).
``provides``
    What the recipe gives to the session, among the same names (e.g.
    ``estimate_auto`` the continuum): the recipes named to the user when a
    prerequisite is missing.
``methods``
    What the recipe did, as a sentence of the methods section of a paper
    (:mod:`astrocook.methods_text`, phase 6-ter PR 5c): a template with the
    parameters as fields (``{kappa}``) and ``{menu}`` where the menu of the
    recipe goes, or a function of the same values that returns the text,
    when the text depends on them.
``references``
    The keys of the works the text of ``methods`` cites, in
    :data:`astrocook.methods_text.REFERENCES`.

Keys of a parameter: ``name``, ``type``, ``default``, ``doc`` (optional,
the name if missing), ``gui_hidden`` (not shown in the dialog), and for the
conversion of its value (:mod:`astrocook.recipes.params`): ``optional``
(None is a value: no value, or automatic), ``keep`` (special strings passed
as they are, e.g. ``'_current_'``), ``min`` and ``max`` (bounds of a
number, included), ``gt`` (a bound excluded: the value must be above it),
``allow_nan`` (NaN is a value; infinities never are), ``choices`` (the
values allowed, for a string: a list in the GUI and for the agents),
``kind`` (what a string names: ``'session'``, an open session,
``'sessions'``, several, comma-separated (phase 6-quater, PR 9: the rule of
both in :mod:`astrocook.recipes._sessions`), or ``'column'``, a column of
the spectrum; phase 6-bis, PR 6c), ``agents`` (False:
not offered to agents, which get the default: the options of debug and the
questions of the GUI, phase 6-bis, PR 6; every recipe is offered, Guido:
the basic ones are blocks of new recipes), ``aliases`` (the old names of a
parameter that was renamed, below).

A renamed parameter keeps its old names as ``aliases`` (phase 6-quater, PR
2, agreed with Guido): the logs saved in the archives, the scripts and the
calls of the MCP tools that use an old name still run, and give the same
result. :meth:`RecipeSpec.canonical` turns the old names into the new one;
it is called where the parameters enter a recipe (the decorator
``typed_recipe``), where the log records them (``resolve_params``), by the
MCP server before it refuses the unknown parameters, and when the log of an
archive is read (:func:`canonical_params`), so that the methods text, the
Log Scripter and the log of the MCP server show the new names only. The
dialogs and the schemas of the MCP tools show the new name only. An alias is
silent (no warning) and permanent: the archives are not changed once saved,
and their logs must replay; an old name given with the new one is an error
(:class:`~astrocook.recipes.errors.InvalidParameter`), even with the same
value. An alias is not in the signature of the method.
"""
from __future__ import annotations

import copy
import logging
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any, Dict, Iterator, Mapping, Optional, Tuple

from astrocook.recipes.errors import InvalidParameter

logger = logging.getLogger(__name__)

#: The categories, in the order of the menus
CATEGORIES: Tuple[str, ...] = ('file', 'edit', 'flux', 'continuum', 'absorbers')

#: What a recipe can require (the key ``requires``)
REQUIREMENTS = frozenset({'z_em', 'resolution', 'continuum'})

_RECIPE_KEYS = {'brief', 'details', 'url', 'params', 'gui_hidden', 'branching',
                'requires', 'provides', 'methods', 'references'}
_PARAM_KEYS = {'name', 'type', 'default', 'doc', 'gui_hidden', 'optional', 'keep',
               'min', 'max', 'gt', 'allow_nan', 'agents', 'choices', 'kind',
               'aliases'}

#: What a string parameter can name (the key ``kind``)
KINDS = frozenset({'session', 'sessions', 'column'})


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
    #: Offered to agents (False: they get the default)
    agents: bool = True
    #: The values allowed (a string parameter), in their spelling
    choices: Tuple[str, ...] = ()
    #: What the value names: ``'session'``, ``'column'``, or "" (nothing)
    kind: str = ""
    #: The old names of the parameter, still accepted (phase 6-quater, PR 2)
    aliases: Tuple[str, ...] = ()
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
    #: What the recipe gives to the session (a subset of :data:`REQUIREMENTS`)
    provides: frozenset = frozenset()
    #: What the recipe did, for the methods text: a template or a function
    #: (:mod:`astrocook.methods_text`); None: a generic sentence
    methods: Any = None
    #: The keys of the works the methods text cites
    references: Tuple[str, ...] = ()

    def param(self, name: str) -> ParamSpec:
        """The parameter called ``name`` (``KeyError`` if there is none)."""
        for p in self.params:
            if p.name == name:
                return p
        raise KeyError(f"{self.category}.{self.name} has no parameter {name!r}")

    def defaults(self) -> Dict[str, Any]:
        """The defaults of all the parameters, by name."""
        return {p.name: p.default for p in self.params}

    def renamed(self) -> Dict[str, str]:
        """The old names of the parameters (their ``aliases``), each with
        its new name."""
        return {old: p.name for p in self.params for old in p.aliases}

    def canonical(self, params: Mapping[str, Any]) -> Dict[str, Any]:
        """``params`` with the old names of the parameters replaced by
        the new ones, in the same order; the other names as they are.
        :class:`InvalidParameter` if a parameter is given by its old and
        its new name."""
        renamed = self.renamed()
        if not renamed or not renamed.keys() & params.keys():
            return dict(params)
        out = {}
        for name, value in params.items():
            new = renamed.get(name)
            if new is None:
                out[name] = value
                continue
            if new in params or new in out:
                raise InvalidParameter(f"{name} is the old name of {new}: give only one "
                                       f"of them")
            logger.debug(f"{self.category}.{self.name}: {name} is the old name of {new}.")
            out[new] = value
        return out


def _param_spec(where: str, entry: Mapping[str, Any]) -> ParamSpec:
    unknown = set(entry) - _PARAM_KEYS
    if unknown:
        raise ValueError(f"{where}: unknown keys {sorted(unknown)}")
    if entry.get('kind', '') not in KINDS | {''}:
        raise ValueError(f"{where}: unknown kind {entry['kind']!r}")
    choices = tuple(entry.get('choices', ()))
    if choices and entry['default'] not in (*choices, ''):  # '' : required
        raise ValueError(f"{where}: the default {entry['default']!r} is not a choice")
    aliases = entry.get('aliases', ())
    if isinstance(aliases, str):
        raise ValueError(f"{where}: the aliases of {entry['name']!r} are a list of names")
    aliases = tuple(aliases)
    for alias in aliases:
        if not isinstance(alias, str) or not alias.isidentifier():
            raise ValueError(f"{where}: the alias {alias!r} is not a name")
    return ParamSpec(name=entry['name'], type=entry['type'],
                     doc=entry.get('doc', entry['name']),
                     gui_hidden=bool(entry.get('gui_hidden', False)),
                     optional=bool(entry.get('optional', False)),
                     keep=tuple(entry.get('keep', ())),
                     min=entry.get('min'), max=entry.get('max'), gt=entry.get('gt'),
                     allow_nan=bool(entry.get('allow_nan', False)),
                     agents=bool(entry.get('agents', True)),
                     choices=choices, kind=str(entry.get('kind', '')),
                     aliases=aliases, _default=copy.deepcopy(entry['default']))


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
    provides = frozenset(entry.get('provides', ()))
    if provides - REQUIREMENTS:
        raise ValueError(f"{where}: unknown provisions {sorted(provides - REQUIREMENTS)}")
    params = tuple(_param_spec(f"{where}", p) for p in entry.get('params', []))
    names = [p.name for p in params]
    if len(set(names)) != len(names):
        raise ValueError(f"{where}: a parameter is declared twice")
    aliases = [a for p in params for a in p.aliases]
    if len(set(aliases)) != len(aliases) or set(aliases) & set(names):
        # An old name means one parameter only
        raise ValueError(f"{where}: an alias is the name or the alias of another "
                         f"parameter")
    return RecipeSpec(category=category, name=name, brief=entry['brief'],
                      details=entry.get('details', ""), url=entry.get('url'),
                      params=params, gui_hidden=bool(entry.get('gui_hidden', False)),
                      branching=bool(entry.get('branching', False)),
                      requires=requires, provides=provides,
                      methods=entry.get('methods'),
                      references=tuple(entry.get('references', ())))


def _build() -> Mapping[str, Mapping[str, RecipeSpec]]:
    from astrocook.recipes.absorbers import ABSORBERS_RECIPES_SCHEMAS
    from astrocook.recipes.continuum import CONTINUUM_RECIPES_SCHEMAS
    from astrocook.recipes.edit import EDIT_RECIPES_SCHEMAS
    from astrocook.recipes.file import FILE_RECIPES_SCHEMAS
    from astrocook.recipes.flux import FLUX_RECIPES_SCHEMAS
    sources: dict[str, Mapping[str, Any]] = {'file': FILE_RECIPES_SCHEMAS, 'edit': EDIT_RECIPES_SCHEMAS,
               'flux': FLUX_RECIPES_SCHEMAS, 'continuum': CONTINUUM_RECIPES_SCHEMAS,
               'absorbers': ABSORBERS_RECIPES_SCHEMAS}
    registry = {}
    seen: dict[str, str] = {}
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


def canonical_params(name: str, params: Mapping[str, Any]) -> Dict[str, Any]:
    """The parameters of a log entry of the recipe ``name``, with the old
    names replaced by the new ones (:meth:`RecipeSpec.canonical`): for the
    log of an archive, when it is read. A recipe that is not in the
    registry, or an entry that gives a parameter by both its names (a log
    written by hand), is kept as it is: the replay says what is wrong."""
    spec = find_recipe(name)
    if spec is None:
        return dict(params)
    try:
        return spec.canonical(params)
    except InvalidParameter:
        return dict(params)
