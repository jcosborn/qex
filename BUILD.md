## QEX application build guide

See [INSTALL.md](INSTALL.md) for installation instructions.

QEX applications can be built using either the standard
`make` command or with `nimble`.
The options and functionality are largely the same, with only a few minor
differences.

## Building executables

| Build Method | Command |
|-|-|
| Make   | make [options] <target> |
| Nimble | nimble make [debug] [FlagsToNim] [Name=Definition] Target [MoreTargets] |



## Building using ``make``

### Quick guide

```
make [options] <target>
```

| Options | Description |
|---------|-------------|
| debug | enable debug build (otherwise release) |
| run | run executable after building |
| verb:N | set build verbosity to integer N |

Target can be the name of a Nim source file (ending with `.nim`)
or a directory (in which case all Nim source in that directory are compiled).

Targets are searched for recursively in the following directories
- `.`
- `qex/src`
- `qex/tests`



## Building using ``nimble``

You can build under the source directory.  Modify the file
`local.nims` to suit your system.

If you wish to build in a separate directory, copy the file
`local.nims` to your build directory, and modify it there.  Copy
`qex.nimble` to your build directory.

Run `nimble tasks` for available tasks, and `nimble help` for
help building your executables.

```
To build nim files:
  nimble make [debug] [FlagsToNim] [Name=Definition] Target [MoreTargets]

`debug' will make debug build.
`Target' can be any file name, optionally with partial directory names.
The produced executables will be under `bin/'.

Examples:
  nimble make debug test0
  nimble make example/testStagProp
```
