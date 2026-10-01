# Score Editor

Tailwind v4 styles and theme tokens live in `input.css` and compile into
`output.css`. The page loads the generated CSS without a Tailwind CDN script.

The editor follows the system's light or dark appearance automatically, including
changes made while the page is open. Both Tailwind utilities and editor syntax
colors use the `prefers-color-scheme` media query.

To rebuild after changing utility classes or the theme:

```sh
cd Score-Editor
npm install
npm run build:css
```

Commit the generated `output.css` along with the source changes so the editor can
also be served directly. Use `npm run watch:css` while editing. The GitHub Pages
workflow rebuilds the CSS before publishing the editor.

If the standalone Tailwind v4 CLI is installed, no npm install is needed:

```sh
tailwindcss -i ./input.css -o ./output.css --minify
```
