from copy import deepcopy
import importlib
from typing import Any, List, Dict

# List of recipes considered 'Branching' or 'Destructive' in Astrocook V2.
# These operations typically involve extraction, combination, or file export,
# requiring the original session state to be preserved and a new entry to be created.
BRANCHING_RECIPES: List[str] = [
    'split',                # Split spectrum by range
    'extract_preset',
    #'part_extract',         # Extracts a sub-part (e.g., specific order/slice)
    'combine',              # Combines multiple sessions
    #'equalize',             # Equalizes multiple sessions
    #'struct_import',        # Imports a structure from another session
    #'save',                 # File output (should always be a new operation)
    #'save_pdf',             # PDF output
    'group_extract',        # Extracts a group into a new structure
    # You will need to add other recipes here as you migrate them (e.g., 'syst_new')
]

def get_recipe_schema(category: str, recipe_name: str) -> Dict[str, Any]:
    """
    Dynamically loads the schemas from the corresponding recipe module.
    :param category: The recipe category (e.g., 'edit', 'fit').
    :param recipe_name: The name of the recipe method (e.g., 'x_convert').
    """
    try:
        # Construct the full path to the module (e.g., 'astrocook.v2.recipes.edit')
        module_path = f"astrocook.v2.recipes.{category}" 
        recipe_module = importlib.import_module(module_path)
        
        # Get the global schemas dictionary from the module (e.g., EDIT_RECIPES_SCHEMAS)
        schemas_dict = getattr(recipe_module, f"{category.upper()}_RECIPES_SCHEMAS")
        
        # Retrieve the specific recipe schema
        schema = schemas_dict[recipe_name]
        
    except (ImportError, AttributeError, KeyError) as e:
        raise NotImplementedError(f"Schema not found for recipe {category}.{recipe_name}: {e}")
        
    return schema

def guarded_deepcopy_v1_state(v1_state_object: Any) -> Any:
    """
    Performs a deepcopy on a V1 state object (like GUILog or Defaults) 
    by temporarily removing the uncopyable '_gui' reference.
    """
    
    # 1. Store and temporarily break the dangerous reference
    original_gui_ref = v1_state_object._gui
    v1_state_object._gui = None

    try:
        # 2. Perform the safe deepcopy
        copied_object = deepcopy(v1_state_object) 
        
    finally:
        # 3. Always restore the reference to the original object
        v1_state_object._gui = original_gui_ref

    # 4. Restore the reference in the copied object as well
    copied_object._gui = original_gui_ref
    
    return copied_object

def is_branching_recipe(recipe_name: str) -> bool:
    """
    Determines if a recipe should create a new session entry 
    (Branching/Destructive operation) or update the current session (Linear operation).
    """
    return recipe_name in BRANCHING_RECIPES