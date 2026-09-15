#!/usr/bin/env python
"""Shared, single-definition detector invocation. Imported by the pilot runner and the tests
so that validation exercises exactly the command the pilot executes."""
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

JAVA_HOME = os.environ.get(
    "FUSIONTAG_JAVA_HOME", r"C:\Program Files\Eclipse Adoptium\jdk-17.0.19.10-hotspot")
P2RANK_HOME = os.path.join(ROOT, "tools", "p2rank_2.5.1")
FPOCKET_WSL = "$HOME/opt/fpocket/bin/fpocket"
WSL_DISTRO = "Ubuntu"

P2RANK_VERSION = "2.5.1"
FPOCKET_VERSION = "source-tag 4.2.3 / commit 4bb0d8447f62fee77e2c3c29f54b5fcaf5e2c066"

# JVM heap taken from the stock prank.bat launcher; not a tuning parameter.
JAVA_OPTS = ["-Xmx2048m"]


def win_to_wsl(path):
    p = os.path.abspath(path).replace("\\", "/")
    return "/mnt/" + p[0].lower() + p[2:] if p[1:3] == ":/" else p


def p2rank_cmd(inp, outdir):
    """P2Rank invoked through java directly (same classpath and main class as prank.bat).

    Calling the `prank` shell script is not portable here: on Windows, `bash` resolves to
    WSL's bash, which cannot execute the Windows-path launcher. Invoking the JVM directly
    removes the shell from the loop entirely.
    """
    cp = os.pathsep.join([os.path.join(P2RANK_HOME, "bin", "p2rank.jar"),
                          os.path.join(P2RANK_HOME, "bin", "lib", "*")])
    return [os.path.join(JAVA_HOME, "bin", "java.exe"), *JAVA_OPTS,
            "-cp", cp, "cz.siret.prank.program.Main",
            "predict", "-f", inp, "-o", outdir, "-threads", "1"]


def p2rank_env():
    return dict(os.environ, JAVA_HOME=JAVA_HOME)


def fpocket_cmd(inp):
    """fpocket through WSL. Returns (argv, env)."""
    argv = ["wsl.exe", "-d", WSL_DISTRO, "--", "bash", "-c",
            f'{FPOCKET_WSL} -f "{win_to_wsl(inp)}"']
    env = dict(os.environ, MSYS2_ARG_CONV_EXCL="*", MSYS_NO_PATHCONV="1")
    return argv, env
