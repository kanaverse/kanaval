include(FetchContent)

FetchContent_Declare(
  millijson
  GIT_REPOSITORY https://github.com/LTLA/millijson
  GIT_TAG 8f2f005fc5b5dc399c7d43d2854c85b42a84ed80
)
FetchContent_MakeAvailable(millijson)
