# philippbroniecki.com

Personal website for Philipp Broniecki — built with [Astro](https://astro.build).

## Quick start

```bash
npm install        # install dependencies
npm run dev        # start dev server at localhost:4321
npm run build      # build for production (output in dist/)
npm run preview    # preview production build locally
```

## Project structure

```
src/
  layouts/
    BaseLayout.astro     # global shell (nav, footer, head, dark mode)
    PageLayout.astro     # content pages with hero heading
    BlogPost.astro       # individual blog post layout
  components/
    Card.astro           # reusable card component (homepage + elsewhere)
    Publication.astro    # academic publication entry with toggle abstract
  pages/
    index.astro          # homepage
    about.astro          # about page (bio, photo, contact links)
    blog/
      index.astro        # blog listing
      [...id].astro      # dynamic blog post route
    academia/
      research.astro     # publications & working papers
      teaching.astro     # teaching history
    rss.xml.ts           # RSS feed
  content/
    blog/                # blog posts (Markdown / MDX)
      welcome.md         # example first post
  styles/
    global.css           # design tokens, reset, utility classes
public/
  images/                # photos, icons
  docs/                  # PDFs (teaching materials, papers)
  data/                  # datasets
  CNAME                  # custom domain for GitHub Pages
```

## Adding a blog post

Create a new `.md` or `.mdx` file in `src/content/blog/`:

```md
---
title: "My new post"
description: "A short description for SEO and the blog index."
pubDate: 2026-04-01
tags: ["data-science", "polling"]
---

Post content in Markdown here.
```

The post will appear at `/blog/<filename>/` and on the blog index automatically.

### Frontmatter fields

| Field         | Required | Description                         |
|---------------|----------|-------------------------------------|
| `title`       | yes      | Post title                          |
| `description` | yes      | Short description (used in SEO/RSS) |
| `pubDate`     | yes      | Publication date (YYYY-MM-DD)       |
| `tags`        | no       | Array of tag strings                |
| `draft`       | no       | Set `true` to hide from listings    |
| `updatedDate` | no       | Last updated date                   |

## Editing homepage content

The homepage is in `src/pages/index.astro`. The main sections are:

- **About** — card linking to the dedicated About page
- **Current Work** — Vantage Insights card, Blog card, and autoMrP card
- **Academia** — cards linking to Research and Teaching pages

The Vantage Insights URL is defined as a variable at the top of the file:
```js
const vantageInsightsUrl = 'https://vantage.philippbroniecki.com';
```
Change this one variable when the domain moves to `insights.philippbroniecki.com` or `insights.vantagedatahouse.com`.

## Editing academic content

- **Research**: `src/pages/academia/research.astro` — uses the `<Publication>` component
- **Teaching**: `src/pages/academia/teaching.astro` — plain structured HTML

## Dark mode

Respects `prefers-color-scheme` by default. Users can toggle manually via the nav button. Preference is saved in `localStorage`.

## URL redirects

Legacy URLs are redirected via `astro.config.mjs`:
- `/research/` → `/academia/research/`
- `/teaching/` → `/academia/teaching/`

## Deployment

The site builds to static HTML (`dist/`). Deploy to GitHub Pages, Netlify, Vercel, or any static host. The `CNAME` file is in `public/` for GitHub Pages custom domain support.

## Tech stack

- **Astro 5** — static site generator
- **Content Collections** — type-safe Markdown/MDX blog
- **@astrojs/sitemap** — auto sitemap generation
- **@astrojs/mdx** — MDX support for rich blog posts
- **@astrojs/rss** — RSS feed
- **CSS custom properties** — theming without a CSS framework
