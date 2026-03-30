# AGENTS.md

Repository guide for agents working on `personal_website`.

## 1. Repository Overview

This repository contains `philippbroniecki.com`, a static personal website built with Astro.

Primary content areas:
- Homepage and About page
- Blog (Markdown/MDX, content collections)
- Academia pages (Research and Teaching)
- RSS feed and sitemap

Key goal when editing: keep content accurate, layout consistent, and builds stable.

## 2. Tech Stack And Runtime

- Astro 5 (`astro`)
- `@astrojs/mdx` for MDX posts/components in content
- `@astrojs/rss` for feed generation
- `@astrojs/sitemap` for sitemap generation
- Static output in `dist/`

Common commands:
- `npm run dev` for local development
- `npm run build` for production build
- `npm run preview` for local preview of built site

## 3. Repository Structure

Important paths:
- `src/pages/` route files
- `src/layouts/` page and post layouts
- `src/components/` reusable UI/components
- `src/content/blog/` blog posts (`.md` and `.mdx`)
- `src/content.config.ts` content collection schema
- `src/styles/global.css` global design tokens and styles
- `public/` static assets copied as-is to output

## 4. Standard Content Workflow

### 4.1 Blog Posts

- Add posts in `src/content/blog/`.
- Use required frontmatter fields from schema:
  - `title`
  - `description`
  - `pubDate`
- Optional fields:
  - `tags`
  - `draft`
  - `updatedDate`

Notes:
- Use `.md` for plain markdown posts.
- Use `.mdx` when importing Astro components into posts.

### 4.2 Site Pages

- Homepage: `src/pages/index.astro`
- About: `src/pages/about.astro`
- Research/Teaching: `src/pages/academia/`

### 4.3 Static Assets

- Place images and files under `public/`.
- Reference with absolute paths like `/images/...`.

## 5. Editing Guardrails

- Prefer minimal, targeted changes.
- Preserve established visual language unless explicitly asked to redesign.
- Do not remove existing functionality while adding new features.
- Keep existing routing and content collection behavior intact.

When possible, verify by running build/dev commands before finalizing.

## 6. Special Section: Iran Interactive Cards In Blog Posts

This section covers embedding the Iran interactive cards with identical style/behavior to the standalone Iran website.

### 6.1 Non-Negotiable Constraint

Do not edit files under:
- `/mnt/storage1/projects/vantage_models/vantage_special_analysis/build_section/`

Treat that location as read-only source artifacts.

### 6.2 Source Of Truth

Canonical Iran standalone artifacts are:
- `/mnt/storage1/projects/vantage_models/vantage_special_analysis/build_section/iran_policy/standalone_website/index.html`
- `/mnt/storage1/projects/vantage_models/vantage_special_analysis/build_section/iran_policy/standalone_website/data/iran_policy_data.json`
- `/mnt/storage1/projects/vantage_models/vantage_special_analysis/build_section/iran_policy/standalone_website/data/iran_policy_crosstabs.json`

Do not rebuild these visuals from scratch in custom Astro components if parity is required.

### 6.3 Integration In This Repository

Current integration points:
- Copied standalone site is served from:
  - `public/iran_policy/standalone_website/`
- Blog embed wrapper component:
  - `src/components/IranCardEmbed.astro`
- Example usage:
  - `src/content/blog/starting-this-blog.mdx`

### 6.4 Refresh Procedure From Upstream Artifacts

From repo root:

```bash
mkdir -p public/iran_policy
rsync -a --delete \
  /mnt/storage1/projects/vantage_models/vantage_special_analysis/build_section/iran_policy/standalone_website/ \
  public/iran_policy/standalone_website/
chmod -R u+rwX,go+rX public/iran_policy/standalone_website
```

### 6.5 Required Embed-Mode Behavior

The copied `public/iran_policy/standalone_website/index.html` must support single-card mode via query params.

Expected param behavior:
- `embed_component=topline_questions`
- `embed_component=demographic_breakdown`
- `embed_component=archetype_breakdown`

Expected aliases supported by the patch:
- `survey_question` -> `topline_questions`
- `demographics` -> `demographic_breakdown`
- `archetypes` -> `archetype_breakdown`

If upstream refresh overwrites local embed-mode patch, re-apply it in the copied `index.html`.

### 6.6 MDX Usage

Use in `.mdx` blog posts:

```mdx
import IranCardEmbed from '../../components/IranCardEmbed.astro';

<IranCardEmbed card="survey_question" />
<IranCardEmbed card="demographics" />
<IranCardEmbed card="archetypes" />
```

Optional:

```mdx
<IranCardEmbed card="demographics" includeToggle={true} />
```

### 6.7 Dark Mode Sync

Blog dark mode and iframe dark mode are independent by default (different localStorage keys).
The integration handles this via two mechanisms:

1. **On iframe load** — `IranCardEmbed.astro` sends a `postMessage` with the current blog theme
   (`{ type: 'vdh_theme', theme: 'dark' | 'light' }`) to every iran card iframe after it loads.
2. **On theme toggle** — a single `MutationObserver` on `document.documentElement` watches
   `data-theme` attribute changes and re-sends the theme to all iframes.
3. **In the iframe** — `public/iran_policy/standalone_website/index.html` has a `message` listener
   (registered only when an `embed_component` URL param is present) that calls `applyTheme()`
   and updates `colorScheme` on `documentElement`.

### 6.8 Transparent Background Fix

In embed mode the standalone site sets `document.body.style.background = 'transparent'`.
Without a matching fix on `<html>`, the browser's canvas color (dark when the OS is in dark mode
and `html` has `color-scheme: light dark`) was visible as a black box around each card.

Fix applied in `applyEmbedComponentMode()`:
```javascript
document.documentElement.style.background = 'transparent';
var resolvedTheme = document.documentElement.getAttribute('data-theme') || 'light';
document.documentElement.style.colorScheme = resolvedTheme === 'dark' ? 'dark' : 'light';
```

If upstream refresh overwrites these lines, re-apply them in the copied `index.html`.

### 6.9 Medium-Width Overflow / Right-Edge Clipping Fix

In embed mode, Demographics and Archetypes can be clipped on the right edge at medium iframe widths
(for example around ~760–840px in blog layouts) if controls/bars force horizontal overflow.

Fixes applied in `public/iran_policy/standalone_website/index.html`:

```css
.layout { grid-template-columns: minmax(0, 1fr); }
.card { min-width: 0; }
```

```css
.embed-component-mode .xt-controls { flex-direction: column; align-items: stretch; }
.embed-component-mode .xt-control-group { min-width: 0; }
.embed-component-mode .xt-select { min-width: 0; width: 100%; max-width: 100%; }
```

```css
@media (max-width: 820px) {
  .embed-component-mode .bar-row  { grid-template-columns: minmax(90px, 8.5rem) minmax(0, 1fr) 66px; }
  .embed-component-mode .mini-bar { grid-template-columns: 40px minmax(0, 1fr) 66px; }
}
```

If upstream refresh overwrites these lines, re-apply them in the copied `index.html`.

### 6.10 Verification Checklist For Iran Card Embeds

- Card style matches Iran standalone website.
- Interactions render correctly in iframe.
- Full standalone chrome is hidden in embedded mode.
- Data files under `public/iran_policy/standalone_website/data/` load correctly.
- No edits were made in `vantage_special_analysis/build_section`.
- Cards have no black box around them in either light or dark mode.
- Cards switch to dark styling when the blog page is switched to dark mode.
- Demographics and Archetypes cards are not cut off on the right side in blog embeds.

## 7. Agent Handoff Expectations

When finishing work, report:
- Files changed
- What behavior changed
- Any verification performed
- Any verification not possible in the environment
