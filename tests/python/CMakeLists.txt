# The 'source' symbolic link simplifies running python tests manually inside
# ${CMAKE_CURRENT_BINARY_DIR}:
# It allows using 'python3 -m pytest [OPTIONS] source/basic_checks.py' instead of
# 'python3 -m pytest [OPTIONS] ../../wsclean/tests/tests/python/basic_checks.py'.
execute_process(
  COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}
          ${CMAKE_CURRENT_BINARY_DIR}/source)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config_vars.py.in
               ${CMAKE_CURRENT_BINARY_DIR}/config_vars.py)

add_custom_target(
  integrationchecks
  COMMAND
    python3 -B -m pytest --exitfirst -v
    --junitxml=${CMAKE_BINARY_DIR}/check_integration.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/basic_checks.py)

add_custom_target(
  facetchecks
  COMMAND
    python3 -B -m pytest --exitfirst -v
    --junitxml=${CMAKE_BINARY_DIR}/check_facets.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/facet_checks.py)

# Test a number of wsclean commands, which can be considered system level tests
add_custom_target(
  longsystemchecks
  COMMAND
    python3 -B -m pytest --exitfirst -v
    --junitxml=${CMAKE_BINARY_DIR}/check_long_system.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vela_deconvolution_checks.py
    ${CMAKE_CURRENT_SOURCE_DIR}/aterms_checks.py
    ${CMAKE_CURRENT_SOURCE_DIR}/long_system_checks.py)
