# t-level

Command line utility for calculating the t-level achieved by running a given number of ECM curves.

## Installation

### Using uv (Recommended)

```bash
uv tool install t-level
```

### Using pip

```bash
pip install t-level
```

## Binary Releases

Standalone binaries for Linux and Windows are available on the [GitHub Releases page](https://github.com/brubsby/t-level/releases).

There are two versions available for each platform:

1.  **Directory Version (Recommended for Speed):**
    *   **Filename:** `t-level-{os}-dir.zip`
    *   **Startup Speed:** Instant (~0.2s).
    *   **Format:** A zip file containing a folder. You must extract it and run the `t-level` executable inside.
    
2.  **Portable Version (Single File):**
    *   **Filename:** `t-level-{os}-portable` (Linux) or `t-level-{os}-portable.exe` (Windows)
    *   **Startup Speed:** Slower (~0.7s delay on every run).
    *   **Format:** A single executable file. Easier to share, but requires unpacking internal data to a temporary directory on every execution.

## Usage

```bash
echo <curve_string>[;<curve_string>][...] | t-level [options]
printf <curve_string>[\n<curve_string>][...] | t-level [options]
t-level [options] < <input_file>
t-level [options] -i <input_file>
t-level [options] -q"<curve_string>[;<curve_string>][...]"
```

`<curve_string>` must fully match the regex:
`\s*(\d+e\d+|\d+)(?:\/(\d+e\d+|\d+))?@(?:B1=)?(\d+e\d+|\d+)(?:,\s*(?:B2=)?(\d+e\d+|\d+))?(?:,\s*(?:(?:param|p)=)?([0-4]))?\s*`

Examples:
- `5208@11e6`
- `5208@11e6,35133391030,1`
- `5208@11e6,35e9,p=3`
- `5208@B1=11e6,B2=35e9,param=0`
- `5208/8192@11e6`  # 5208 finished curves, 2984 stage-1 curves

Multiple curve strings must be delimited by semicolons or newlines.

## Options

See `t-level --help` for more details.

