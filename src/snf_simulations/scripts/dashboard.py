"""Launch script for the SNF antineutrino spectrum dashboard."""

from snf_simulations.dashboard import app


def main() -> None:
    """Run the Shiny dashboard."""
    app.run()


if __name__ == "__main__":
    main()
