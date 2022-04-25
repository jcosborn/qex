## Building using ``make``

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
