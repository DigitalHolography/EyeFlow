"""Frozen launcher, including the build-time runtime verification command."""
import sys

from launcher import main


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--packaging-check":
        from frozen_check import main as check_main
        raise SystemExit(check_main(sys.argv[2:]))
    main()
