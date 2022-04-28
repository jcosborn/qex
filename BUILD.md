## QEX application build guide

See [INSTALL.md](INSTALL.md) for installation instructions.

QEX applications can be built using either the standard
`make` command or with `nimble`.
The options and functionality are largely the same, with only a few minor
differences.

### Building executables

| Build Method | Command |
| ---- | --- |
| Make   | make [build option \| Nim option]... [path]... |
| Nimble | nimble make [build option \| Nim option]... [path]... |

Details on the build options are [below](#build-options).

The `path` is the name of the target to be compiled.
It can include wildcards and parts of the path.
Targets are searched for matches within the current directory
and the nim source (`src`) and tests directories.

You can see a list of matching targets with the `targets` command, e.g.
```
make targets <path>
```
or
```
nimble targets <path>
```


## Full command syntax

| Build Method | Command |
| ---- | --- |
| Make   | make [COMMAND] [build option \| Nim option]... [path]... |
| Nimble | nimble COMMAND [build option \| Nim option]... [path]... |

The `command` is required for nimble, but optional for make (it will default
to the `make` command below).

### Commands

| Command | Description |
| --- | --- |
| help    | Show help message |
| show    | Show Nim compile flags |
| targets | Show available build targets. Targets <name> will search for targets matching <name> (can include standard shell wildcards) |
| clean   | Remove contents of nimcache directory |
| tests   | Build tests and create `testscript.sh` test runner |
| make    | Search for each [path]... as described below, compile, link, and put executables in `bin` |

When using the `make` build method the `make` command is default and can be skipped.


### Build options

| Option | Description |
| --- | --- |
|  cc    | compile in C mode |
|  cpp   | compile in C++ mode |
|  debug | set debug build |
|  run   | run executable after building |
|  verb  | set build verbosity to N (verb:N), N in 0,1,2,3 |


### Nim options

| Option | Description |
| --- | --- |
| -<option>  | Passes `-<option>` to Nim compiler |
| :-<option> | Passes `-<option>` to Nim compiler (avoids issues with make/nimble trying to parse it). |
| :foo       | Sets Nim define `foo` (equivalent to `-d:foo`). |
| :foo=bar   | Sets Nim define `foo` to value `bar` (equivalent to `-d:foo=bar`). |

### Path

| PATH | Description |
| --- | --- |
| foo.nim  | Search for file matching `*foo.nim` in source paths (including subdirectories, but not following links) |
|  foo     | Search for both `*foo.nim` and `*foo`, if a directory matches compile all `*.nim` in it |

Note:  only one match is allowed,
       specify part of path to resolve ambiguity

source paths: `.`, `qex/src`, `qex/tests`


### Examples:

```
  make debug test0
  make example/testStagProp
```

```
  nimble make debug test0
  nimble make example/testStagProp
```
