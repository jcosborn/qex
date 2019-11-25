type
  StdString* {.importcpp: "std::string", header: "<string>".} = object

proc newStdString*(): StdString {.
  constructor, importCpp: "std::string()", header: "<string>".}

proc newStdString*(s: cstring): StdString {.
  constructor, importCpp: "std::string(@)", header: "<string>".}

proc size*(this: var StdString): csize {.
  importCpp: "size", header: "<string>".}

proc append*(this: var StdString, str: StdString): var StdString {.
  importcpp: "#.append(#)", header: "<string>".}

proc cstr*(x: StdString): cstring {.
  importcpp: "#.c_str()", header: "<string>".}

proc `$`*(x: StdString): string =
  $(cstr(x))

# Streams

type
  IStream* {.header: "<istream>", importcpp: "std::istream".} = object
  IStringStream* {.header: "<sstream>", importcpp: "std::istringstream".} = object

proc newIStringStream*(foo: cstring): IStringStream {.
  constructor, importCpp: "std::istringstream(@)", header: "<sstream>"}


when isMainModule:
  var
    s1 = newStdString("test1")
    s2 = newStdString("test2")

  echo s1.append s2

  var myString = newStdString()
  echo myString.size()
  echo myString

  let s = newStdString("test")
  echo s
